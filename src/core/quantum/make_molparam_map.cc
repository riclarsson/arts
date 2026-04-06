#include <cstdlib>
#include <filesystem>
#include <ostream>

#include "hitran_species_info.h"

namespace {
std::vector<HitranSpeciesInfo> read(std::filesystem::path path) {
  std::vector<HitranSpeciesInfo> out;

  for (const auto& entry : std::filesystem::directory_iterator(path)) {
    if (entry.is_regular_file() && entry.path().extension() == ".xml") {
      xml_read_from_file(entry.path().string(), out.emplace_back());
    }
  }

  return out;
}

std::map<Index, std::map<char, std::pair<Index, Numeric>>> from(
    const std::vector<HitranSpeciesInfo>& data) {
  std::map<Index, std::map<char, std::pair<Index, Numeric>>> out;

  for (const auto& x : data) {
    out[x.hitind][x.hitchar] = {Species::find_species_index(x.spec), x.ratio};
  }

  return out;
}
}  // namespace

void write_header(
    std::ostream& os,
    const std::map<Index, std::map<char, std::pair<Index, Numeric>>>& data) {
  std::println(os, R"(#pragma once

#include <map>
#include <utility>

#include <isotopologues.h>

namespace Hitran {{
std::map<Index, std::map<char, std::pair<Index, Numeric>>> molparam_map{{)");

  for (const auto& [mol, iso_map] : data) {
    std::print(os, "    {{{}, {{", mol);
    for (const auto& [isochar, pair] : iso_map) {
      std::print(
          os, "{{'{}', {{{}, {}}}}}, ", isochar, pair.first, pair.second);
    }
    std::println(os, "}}}},");
  }

  std::println(os, R"(}};
}}  // namespace Hitran)");
}

int main() try {
  std::filesystem::path path = "/Users/richard/Work/arts-cat-data/";
  auto data                  = read(path / "hitran/");

  std::ofstream h("hitran_species_map.h", std::ios::out);

  write_header(h, from(std::move(data)));

  return EXIT_SUCCESS;
} catch (const std::exception& e) {
  std::println(stderr, "Error in main: {}", e.what());
  return EXIT_FAILURE;
}