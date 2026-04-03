#include <filesystem>
#include <print>
#include <set>

#include "species_enum_info.h"

namespace {
std::set<SpeciesEnumInfo> read_split_species(
    const std::filesystem::path& path) {
  std::set<SpeciesEnumInfo> species_set{};

  for (const auto& entry : std::filesystem::directory_iterator(path)) {
    if (entry.is_regular_file() && entry.path().extension() == ".xml") {
      SpeciesEnumInfo info;
      xml_read_from_file(entry.path().string(), info);

      if (species_set.count(info) > 0) {
        throw std::runtime_error(std::format(
            "Duplicate species enum info found: {}", info.shortname));
      }

      species_set.insert(info);
    }
  }

  return species_set;
}

void write_header(std::ostream& os) {
  std::println(os,
               R"(#pragma once

#include <species_enum_info.h>

#include <vector>

const std::vector<SpeciesEnumInfo>& get_species_enum_info();
)");
}

void write_source(std::ostream& os,
                  const std::set<SpeciesEnumInfo>& species_set) {
  std::println(os,
               R"(#include "auto_species_enum_info.h"

namespace {{
std::vector<SpeciesEnumInfo> get_species_enum_info_local() {{
  std::vector<SpeciesEnumInfo> info{{}};
  info.reserve({});
)",
               species_set.size());

  for (const auto& info : species_set) {
    std::println(os,
                 R"(  info.emplace_back({}, "{}", "{}");)",
                 info.enum_value,
                 info.shortname,
                 info.longname);
  }

  std::println(os, R"(
  return info;
}}
}}  // namespace

const std::vector<SpeciesEnumInfo>& get_species_enum_info() {{
  const static auto& info = get_species_enum_info_local();
  return info;
}})");
}
}  // namespace

int main(int argc, char** argv) try {
  if (argc != 2) {
    std::println(stderr, "Usage: {} <species>", argv[0]);
    return EXIT_FAILURE;
  }

  const std::filesystem::path path = argv[1];

  auto h  = std::fstream("auto_species_enum_info.h", std::ios::out);
  auto cc = std::fstream("auto_species_enum_info.cc", std::ios::out);

  write_source(cc, read_split_species(path / "species/"));
  write_header(h);

} catch (const std::exception& e) {
  std::println(stderr, "Error: {}", e.what());
  return EXIT_FAILURE;
} catch (...) {
  std::println(stderr, "An unknown error occurred.");
  return EXIT_FAILURE;
}
