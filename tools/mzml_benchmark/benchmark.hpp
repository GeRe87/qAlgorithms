#ifndef QALGORITHMS_MZML_BENCHMARK_HPP
#define QALGORITHMS_MZML_BENCHMARK_HPP

#include <cerrno>
#include <filesystem>
#include <string>
#include <system_error>
#include <vector>
#ifdef _WIN32
#include <process.h>
#else
#include <spawn.h>
#include <sys/wait.h>
extern char **environ;
#endif

namespace qAlgorithms::benchmark
{
    // The arguments are identical to the CLI. Results are written to the
    // requested paths, and the return value is the process exit code.
    inline int run(const std::filesystem::path &python, const std::filesystem::path &script,
                   const std::vector<std::filesystem::path> &arguments)
    {
        std::vector<std::filesystem::path::string_type> command{python.native(), script.native()};
        for (const auto &argument : arguments)
            command.push_back(argument.native());
#ifdef _WIN32
        // The Windows spawn API concatenates its arguments before passing
        // them to the child. Quote each entry for the Windows argv parser.
        for (auto &argument : command)
        {
            std::wstring escaped(1, L'"');
            size_t backslashes = 0;
            for (const auto c : argument)
            {
                if (c == L'\\')
                {
                    ++backslashes;
                    continue;
                }
                escaped.append(c == L'"' ? 2 * backslashes + 1 : backslashes, L'\\');
                escaped += c;
                backslashes = 0;
            }
            escaped.append(2 * backslashes, L'\\');
            argument = escaped + L'"';
        }
        std::vector<const wchar_t *> argv;
        for (const auto &argument : command)
            argv.push_back(argument.c_str());
        argv.push_back(nullptr);
        const auto status = _wspawnvp(_P_WAIT, python.c_str(), argv.data());
        if (status == -1)
            throw std::system_error(errno, std::generic_category(), "Could not start Python");
        return static_cast<int>(status);
#else
        std::vector<char *> argv;
        for (auto &argument : command)
            argv.push_back(argument.data());
        argv.push_back(nullptr);
        pid_t child;
        const int error = posix_spawnp(&child, python.c_str(), nullptr, nullptr, argv.data(), environ);
        if (error != 0)
            throw std::system_error(error, std::generic_category(), "Could not start Python");
        int status;
        while (waitpid(child, &status, 0) == -1)
            if (errno != EINTR)
                throw std::system_error(errno, std::generic_category(), "Could not wait for Python");
        return WIFEXITED(status) ? WEXITSTATUS(status) : 1;
#endif
    }
}
#endif
