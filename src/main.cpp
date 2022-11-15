// This file will be generated automatically when you run the CMake configuration step.
// It creates a namespace called `myproject`.
// You can modify the source template at `configured_files/config.hpp.in`.
#include <internal_use_only/config.hpp>

#include <fmt/format.h>


// NOLINTNEXTLINE(bugprone-exception-escape)
int main()
{
  fmt::print("{}\n", myproject::cmake::project_version);
  return EXIT_SUCCESS;
}
