#include "xml_subsurface.h"

void xml_io_stream<SubsurfaceField>::write(std::ostream& os,
                                           const SubsurfaceField& x,
                                           bofstream* pbofs,
                                           std::string_view name) {
  std::println(os, R"(<{0} name="{1}">)", type_name, name);

  xml_write_to_stream(os, x.bottom_depth, pbofs, "Depth"sv);
  xml_write_to_stream(os, x.map_data, pbofs);

  std::println(os, R"(</{0}>)", type_name);
}

void xml_io_stream<SubsurfaceField>::read(std::istream& is,
                                          SubsurfaceField& x,
                                          bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, x.bottom_depth, pbifs);
  xml_read_from_stream(is, x.map_data, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

void xml_io_stream<SubsurfaceData>::read(std::istream& is,
                                         SubsurfaceData& v,
                                         bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, v.data, pbifs);
  xml_read_from_stream(is, v.alt_upp, pbifs);
  xml_read_from_stream(is, v.alt_low, pbifs);
  xml_read_from_stream(is, v.lat_upp, pbifs);
  xml_read_from_stream(is, v.lat_low, pbifs);
  xml_read_from_stream(is, v.lon_upp, pbifs);
  xml_read_from_stream(is, v.lon_low, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

void xml_io_stream<SubsurfaceData>::write(std::ostream& os,
                                          const SubsurfaceData& v,
                                          bofstream* pbofs,
                                          std::string_view name) {
  std::println(os, R"(<{0} name="{1}">)", type_name, name);

  xml_write_to_stream(os, v.data, pbofs, "Data"sv);
  xml_write_to_stream(os, v.alt_upp, pbofs, "alt_upp"sv);
  xml_write_to_stream(os, v.alt_low, pbofs, "alt_low"sv);
  xml_write_to_stream(os, v.lat_upp, pbofs, "lat_upp"sv);
  xml_write_to_stream(os, v.lat_low, pbofs, "lat_low"sv);
  xml_write_to_stream(os, v.lon_upp, pbofs, "lon_upp"sv);
  xml_write_to_stream(os, v.lon_low, pbofs, "lon_low"sv);

  std::println(os, R"(</{0}>)", type_name);
}

void xml_io_stream<SubsurfacePoint>::read(std::istream& is,
                                          SubsurfacePoint& v,
                                          bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, v.temperature, pbifs);
  xml_read_from_stream(is, v.density, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

void xml_io_stream<SubsurfacePoint>::write(std::ostream& os,
                                           const SubsurfacePoint& v,
                                           bofstream* pbofs,
                                           std::string_view name) {
  std::println(os, R"(<{0} name="{1}">)", type_name, name);

  xml_write_to_stream(os, v.temperature, pbofs, "temperature"sv);
  xml_write_to_stream(os, v.density, pbofs, "density"sv);

  std::println(os, R"(</{0}>)", type_name);
}
