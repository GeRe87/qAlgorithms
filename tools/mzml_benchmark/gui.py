"""Small desktop front end for the benchmark commands."""

from pathlib import Path
import queue
import threading
import tkinter as tk
from tkinter import filedialog, ttk

from benchmark import build as build_benchmark, default_profiler, save, score


def launch():
    window = tk.Tk()
    window.title("qAlgorithms · mzML-Benchmark")
    window.columnconfigure(0, weight=1)
    tabs = ttk.Notebook(window)
    tabs.grid(row=0, column=0, padx=12, pady=12, sticky="nsew")
    status = tk.StringVar(value="Bereit")
    ttk.Label(window, textvariable=status, wraplength=620).grid(row=1, column=0, padx=12, pady=(0, 12), sticky="w")
    jobs = queue.Queue()
    buttons = []

    def submit(task):
        for button in buttons:
            button.state(["disabled"])
        status.set("Verarbeitung läuft …")

        def run():
            try:
                jobs.put(task())
            except Exception as error:
                jobs.put(f"Fehler: {error}")

        threading.Thread(target=run, daemon=True).start()

    def poll():
        try:
            status.set(jobs.get_nowait())
        except queue.Empty:
            pass
        else:
            for button in buttons:
                button.state(["!disabled"])
        window.after(100, poll)

    def panel(title):
        frame = ttk.Frame(tabs, padding=14)
        frame.columnconfigure(1, weight=1)
        tabs.add(frame, text=title)
        return frame

    def field(frame, row, label, default="", picker=None):
        value = tk.StringVar(value=default)
        ttk.Label(frame, text=label).grid(row=row, column=0, sticky="w", padx=(0, 10), pady=5)
        ttk.Entry(frame, textvariable=value, width=54).grid(row=row, column=1, sticky="ew", pady=5)
        if picker:
            def select():
                chosen = picker()
                if chosen:
                    value.set(chosen)
            ttk.Button(frame, text="Auswählen", command=select).grid(row=row, column=2, padx=(8, 0))
        return value

    def action(frame, row, title, prepare):
        def click():
            try:
                submit(prepare())
            except (ValueError, OSError) as error:
                status.set(f"Fehler: {error}")
        button = ttk.Button(frame, text=title, command=click)
        button.grid(row=row, column=1, sticky="w", pady=(10, 0))
        buttons.append(button)

    build = panel("Benchmark erstellen")
    inputs = ()
    picked = tk.StringVar(value="Keine Dateien ausgewählt")

    def choose_inputs():
        nonlocal inputs
        selected = filedialog.askopenfilenames(filetypes=[("mzML", "*.mzML"), ("Alle Dateien", "*")])
        if selected:
            inputs = selected
            picked.set(f"{len(inputs)} Datei(en) ausgewählt")

    ttk.Label(build, textvariable=picked).grid(row=0, column=1, sticky="w")
    ttk.Button(build, text="mzML-Dateien …", command=choose_inputs).grid(row=0, column=2)
    output = field(build, 1, "Ausgabeordner", picker=filedialog.askdirectory)
    profiler = field(build, 2, "qAlgorithms-Profiler", default_profiler(), filedialog.askopenfilename)
    seed = field(build, 3, "Seed", "1")
    peaks = field(build, 4, "Features pro Run", "12")
    runs = field(build, 5, "Synthetische Runs", "1")
    kernel = tk.StringVar(value="quadratic")
    ttk.Label(build, text="Kalibrierungskernel").grid(row=6, column=0, sticky="w")
    ttk.Combobox(build, textvariable=kernel, values=("quadratic", "hyperbolic"),
                 state="readonly").grid(row=6, column=1, sticky="w", pady=5)
    replicate_inputs = tk.BooleanVar()
    ttk.Checkbutton(build, text="Eingabedateien sind Replikate derselben Probe",
                    variable=replicate_inputs).grid(row=7, column=1, sticky="w")

    def prepare_build():
        if not inputs or not output.get():
            raise ValueError("Wähle mindestens eine mzML-Datei und einen Ausgabeordner.")
        sources, target, native = inputs, Path(output.get()), profiler.get()
        used_seed, count, repeat = int(seed.get()), int(peaks.get()), int(runs.get())
        is_replicates = replicate_inputs.get()
        model = kernel.get()
        if count < 0 or repeat < 1:
            raise ValueError("Wähle eine nichtnegative Feature-Zahl und mindestens einen synthetischen Run.")

        def task():
            build_benchmark(sources, target, profiler=native, kernel=model, replicate_inputs=is_replicates,
                            seed=used_seed, peaks=count, replicates=repeat)
            return f"Benchmark erstellt. {target / 'synthetic'}"
        return task

    action(build, 8, "Kalibrieren und erzeugen", prepare_build)
    assess = panel("Feature-Liste bewerten")
    truth = field(assess, 0, "Ground Truth (JSON)", picker=filedialog.askopenfilename)
    found = field(assess, 1, "Feature-Liste (CSV)", picker=filedialog.askopenfilename)
    ppm = field(assess, 2, "m/z-Toleranz (ppm)", "5")
    seconds = field(assess, 3, "RT-Toleranz (s)", "5")
    score_output = field(assess, 4, "Ergebnisdatei (optional)",
                         picker=lambda: filedialog.asksaveasfilename(defaultextension=".json"))

    def prepare_score():
        expected, detected = truth.get(), found.get()
        if not expected or not detected:
            raise ValueError("Wähle Ground Truth und Feature-Liste.")
        mass_tol, rt_tol = float(ppm.get()), float(seconds.get())
        target = score_output.get()

        def task():
            result = score(expected, detected, mass_tol, rt_tol)
            if target:
                save(Path(target), result)
            return f"TP {result['tp']} · FP {result['fp']} · FN {result['fn']}"
        return task

    action(assess, 5, "Bewerten", prepare_score)
    poll()
    window.mainloop()
