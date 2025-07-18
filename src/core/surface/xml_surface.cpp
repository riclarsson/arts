#include "xml_surface.h"

void xml_io_stream<SurfaceField>::write(std::ostream& os,
                                        const SurfaceField& x,
                                        bofstream* pbofs,
                                        std::string_view name) {
  std::println(os, R"(<{0} name="{1}">)", type_name, name);

  xml_write_to_stream(os, x.ellipsoid, pbofs, "Ellipsoid"sv);
  xml_write_to_stream(os, x.map_data, pbofs);

  std::println(os, R"(</{0}>)", type_name);
}

void xml_io_stream<SurfaceField>::read(std::istream& is,
                                       SurfaceField& x,
                                       bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, x.ellipsoid, pbifs);
  xml_read_from_stream(is, x.map_data, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

void xml_io_stream<SurfaceData>::write(std::ostream& os,
                                       const SurfaceData& x,
                                       bofstream* pbofs,
                                       std::string_view name) {
  std::println(os, R"(<{0} name="{1}">)", type_name, name);

  xml_write_to_stream(os, x.data, pbofs);
  xml_write_to_stream(os, x.lat_upp, pbofs);
  xml_write_to_stream(os, x.lat_low, pbofs);
  xml_write_to_stream(os, x.lon_upp, pbofs);
  xml_write_to_stream(os, x.lon_low, pbofs);

  std::println(os, R"(</{0}>)", type_name);
}

void xml_io_stream<SurfaceData>::read(std::istream& is,
                                      SurfaceData& x,
                                      bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, x.data, pbifs);
  xml_read_from_stream(is, x.lat_upp, pbifs);
  xml_read_from_stream(is, x.lat_low, pbifs);
  xml_read_from_stream(is, x.lon_upp, pbifs);
  xml_read_from_stream(is, x.lon_low, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

void xml_io_stream<SurfacePoint>::write(std::ostream& os,
                                        const SurfacePoint& x,
                                        bofstream* pbofs,
                                        std::string_view name) {
  std::println(os, R"(<{0} name="{1}">)", type_name, name);

  xml_write_to_stream(os, x.elevation, pbofs);
  xml_write_to_stream(os, x.temperature, pbofs);
  xml_write_to_stream(os, x.normal, pbofs);
  xml_write_to_stream(os, x.prop, pbofs);

  std::println(os, R"(</{0}>)", type_name);
}

void xml_io_stream<SurfacePoint>::read(std::istream& is,
                                       SurfacePoint& x,
                                       bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, x.elevation, pbifs);
  xml_read_from_stream(is, x.temperature, pbifs);
  xml_read_from_stream(is, x.normal, pbifs);
  xml_read_from_stream(is, x.prop, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

void xml_io_stream<SurfacePropertyTag>::write(std::ostream& os,
                                              const SurfacePropertyTag& x,
                                              bofstream* pbofs,
                                              std::string_view name) {
  std::println(os, R"(<{0} name="{1}">)", type_name, name);

  xml_write_to_stream(os, x.name, pbofs);

  std::println(os, R"(</{0}>)", type_name);
}

void xml_io_stream<SurfacePropertyTag>::read(std::istream& is,
                                             SurfacePropertyTag& x,
                                             bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, x.name, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}
