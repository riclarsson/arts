////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////
/*!
  \file   xml_io.cc
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2002-05-10

  \brief This file contains basic functions to handle XML data files.

*/

#include "xml_io_base.h"

#include <bifstream.h>
#include <bofstream.h>
#include <double_imanip.h>
#include <file.h>

#include <format>
#include <iomanip>
#include <iterator>
#include <stdexcept>
#include <string_view>

constexpr Index ARTS_XML_VERSION = 3;

namespace {
static inline std::string quotation_mark_replacement{"”"};
static inline std::string quotation_mark_original{"\""};
}  // namespace

////////////////////////////////////////////////////////////////////////////
//   XMLTag implementation
////////////////////////////////////////////////////////////////////////////

//! Check tag name
/*!
  Checks whether the name of the tag is correct. Throws runtime
  error otherwise.

  \param expected_name Expected tag name
*/
void XMLTag::check_name(const std::string_view& expected_name) {
  if (name != expected_name)
    xml_parse_error(
        std::format("Tag <{}> expected but <{}> found.", expected_name, name));
}
void XMLTag::check_end_name(const std::string_view& expected_name) {
  if (name.empty()) {
    xml_parse_error(std::format(
        "Tag </{}> expected but nothing found.", expected_name, name));
  }

  if (name.front() != '/' or
      std::string_view{name.begin() + 1, name.end()} != expected_name)
    xml_parse_error(
        std::format("Tag </{}> expected but <{}> found.", expected_name, name));
}

/*! Adds a String attribute to tag

All " are replaced by ” to work in the XML tag

  \param aname Attribute name
  \param value Attribute value
*/
void XMLTag::add_attribute(const std::string_view& aname,
                           const std::string_view& value_) {
  String value{value_};

  String name{aname};

  auto pos = value.find(quotation_mark_original);
  while (pos not_eq std::string::npos) {
    value.replace(
        pos, quotation_mark_original.length(), quotation_mark_replacement);
    pos = value.find(quotation_mark_original);
  }

  attribs[name] = value;
}

//! Adds an Index attribute to tag
/*!

  \param aname Attribute name
  \param value Attribute value
*/
void XMLTag::add_attribute(const std::string_view& aname, const Index& value) {
  add_attribute(aname, std::format("{}", value));
}

void XMLTag::add_attribute(const std::string_view& aname, const Size& value) {
  add_attribute(aname, std::format("{}", value));
}

void XMLTag::add_attribute(const std::string_view& aname,
                           const Numeric& value) {
  add_attribute(aname, std::format("{}", value));
}

//! Checks whether attribute has the expected value
/*!

  If the attribute has another value or is unknown an exception is
  thrown.

  \param aname Attribute name
  \param value Expected value
*/
void XMLTag::check_attribute(const std::string_view& aname,
                             const std::string_view& value) {
  String actual_value;

  get_attribute_value(aname, actual_value);

  if (actual_value == "*not found*") {
    xml_parse_error(std::format("Required attribute {} does not exist", aname));
  } else if (actual_value != value) {
    xml_parse_error(
        std::format("Attribute {} has value \"{}\" but \"{}\" was expected.",
                    aname,
                    actual_value,
                    value));
  }
}

void XMLTag::check_attribute(const std::string_view& aname, const Size& value) {
  Size actual_value;

  get_attribute_value(aname, actual_value);

  if (actual_value != value) {
    xml_parse_error(
        std::format("Attribute {} has value \"{}\" but \"{}\" was expected.",
                    aname,
                    actual_value,
                    value));
  }
}

void XMLTag::check_attribute(const std::string_view& aname,
                             const Index& value) {
  Index actual_value;

  get_attribute_value(aname, actual_value);

  if (actual_value != value) {
    xml_parse_error(
        std::format("Attribute {} has value \"{}\" but \"{}\" was expected.",
                    aname,
                    actual_value,
                    value));
  }
}

bool XMLTag::has_attribute(const std::string_view& aname) const {
  return attribs.contains(String{aname});
}

/*! Returns value of attribute as String

  Searches for the matching attribute and returns it value. If no
  attribute with the given name exists, return value is set to
  *not found*.

  Replaces all ” with " to counter XML-tag problems

  \param aname Attribute name
  \param value Return value
*/
void XMLTag::get_attribute_value(const std::string_view& aname,
                                 String& value,
                                 std::string_view def) try {
  value = attribs.at(String{aname});

  auto pos = value.find(quotation_mark_replacement);
  while (pos not_eq std::string::npos) {
    value.replace(
        pos, quotation_mark_replacement.length(), quotation_mark_original);
    pos = value.find(quotation_mark_replacement);
  }
} catch (const std::out_of_range&) {
  // If the attribute does not exist, set value to *not found*
  if (def.empty())
    value = std::format("*{} not found*", aname);
  else
    value = String{def};
  // and do not throw an exception
}

//! Returns value of attribute as type Index
/*!
  Searches for the matching attribute and returns it value. If no
  attribute with the given name exists, return value is set to
  *not found*.

  \param aname Attribute name
  \param value Return value
*/
void XMLTag::get_attribute_value(const std::string_view& aname, Index& value) {
  String attribute_value;
  std::istringstream strstr("");

  get_attribute_value(aname, attribute_value);
  strstr.str(attribute_value);
  strstr >> value;
  if (strstr.fail()) {
    xml_parse_error(std::format(
        "Error while parsing value of {0} from <{1}>", aname, name));
  }
}

void XMLTag::get_attribute_value(const std::string_view& aname, Size& value) {
  String attribute_value;
  std::istringstream strstr("");

  get_attribute_value(aname, attribute_value);
  strstr.str(attribute_value);
  strstr >> value;
  if (strstr.fail()) {
    xml_parse_error(std::format(
        "Error while parsing value of {0} from <{1}>", aname, name));
  }
}

void XMLTag::get_attribute_value(const std::string_view& aname,
                                 Numeric& value) {
  String attribute_value;
  std::istringstream strstr("");

  get_attribute_value(aname, attribute_value);
  strstr.str(attribute_value);
  strstr >> double_imanip() >> value;
  if (strstr.fail()) {
    xml_parse_error(std::format(
        "Error while parsing value of {0} from <{1}>", aname, name));
  }
}

//! Reads next XML tag
/*!
  Reads the name and attributes of the next XML tag from stream.

  \param is Input stream
*/
void XMLTag::read_from_stream(std::istream& is) {
  String token, attr_name, attr_value;
  std::stringbuf tag;
  std::istringstream sstr("");

  char ch = 0;

  attribs.clear();

  while (is.good() && isspace(is.peek())) {
    is.get();
  }

  if (not is.good()) xml_parse_error("Unknown error while reading from stream");

  is >> ch;

  if (ch != '<') {
    is >> token;

    if (ch != '\0') token = ch + token;

    xml_parse_error(std::format("'<' expected but '{}' found.", token));
  }

  is.get(tag, '>');

  // Hit EOF while looking for '>'
  if (is.bad() || is.eof()) {
    xml_parse_error("Unexpected end of file while looking for '>'");
  }

  if (is.get() != '>') {
    xml_parse_error("Closing > not found in tag: " + tag.str());
  }

  sstr.str(tag.str() + '>');

  sstr >> name;

  if (name[name.length() - 1] == '>') {
    // Because closin > was found, the tag for sure has no
    // attributes, set token to ">" to skip reading of attributes
    name.erase(name.length() - 1, 1);
    token = ">";
  } else {
    // Tag may have attributes, so read next token
    sstr >> token;
  }

  //extract attributes
  while (token != ">") {
    String::size_type pos;

    pos = token.find('=', 0);
    if (pos == String::npos) {
      xml_parse_error("Syntax error in tag: " + tag.str());
    }

    attr_name = token.substr(0, pos);
    token.erase(0, pos + 1);

    if (token[0] != '\"') {
      xml_parse_error("Missing \" in tag: " + tag.str());
    }

    while ((pos = token.find('\"', 1)) == (String::size_type)String::npos &&
           token != ">") {
      String ntoken;
      sstr >> ntoken;
      if (!ntoken.length()) break;
      token += " " + ntoken;
    }

    if (pos == (String::size_type)String::npos) {
      xml_parse_error("Missing \" in tag: " + sstr.str());
    }

    if (pos == 1)
      attr_value = "";
    else
      attr_value = token.substr(1, pos - 1);

    attribs[attr_name] = attr_value;

    if (token[token.length() - 1] == '>') {
      token = ">";
    } else {
      sstr >> token;
    }
  }

  // Skip comments
  if (name == "comment") {
    is.get(tag, '<');

    // Hit EOF while looking for '<'
    if (is.bad() || is.eof()) {
      xml_parse_error(
          "Unexpected end of file while looking for "
          "comment tag");
    }

    read_from_stream(is);
    check_name("/comment"sv);
    read_from_stream(is);
  }
}

//! Write XML tag
/*!
  Puts the tag together and writes it to stream.

  \param os Output stream
*/
void XMLTag::write_to_stream(std::ostream& os) {
  std::print(os, "<{}", name);

  for (auto& [an, av] : attribs) std::print(os, R"( {}="{}")", an, av);

  std::print(os, ">");
}

////////////////////////////////////////////////////////////////////////////
//   Functions to open and read XML files
////////////////////////////////////////////////////////////////////////////

//! Open file for XML output
/*!
  This function opens an XML file for writing.

  \param file Output filestream
  \param name Filename
*/
void xml_open_output_file(std::ofstream& file, const String& name) {
  // Tell the stream that it should throw exceptions.
  // Badbit means that the entire stream is corrupted, failbit means
  // that the last operation has failed, but the stream is still
  // valid. We don't want either to happen!
  // FIXME: This does not yet work in  egcs-2.91.66, try again later.
  file.exceptions(std::ios::badbit | std::ios::failbit);

  // c_str explicitly converts to c String.
  try {
    file.open(name.c_str());
  } catch (const std::exception&) {
    std::ostringstream os;
    os << "Cannot open output file: " << name << '\n'
       << "Maybe you don't have write access "
       << "to the directory or the file?";
    throw std::runtime_error(os.str());
  }

  // See if the file is ok.
  // FIXME: This should not be necessary anymore in the future, when
  // g++ stream exceptions work properly. (In that case we would not
  // get here if there really was a problem, because of the exception
  // thrown by open().)
  if (!file) {
    std::ostringstream os;
    os << "Cannot open output file: " << name << '\n'
       << "Maybe you don't have write access "
       << "to the directory or the file?";
    throw std::runtime_error(os.str());
  }
}

#ifdef ENABLE_ZLIB

//! Open file for zipped XML output
/*!
  This function opens a zipped XML file for writing.

  \param file Output filestream
  \param name Filename
*/
void xml_open_output_file(ogzstream& file, const String& name) {
  // Tell the stream that it should throw exceptions.
  // Badbit means that the entire stream is corrupted, failbit means
  // that the last operation has failed, but the stream is still
  // valid. We don't want either to happen!
  // FIXME: This does not yet work in  egcs-2.91.66, try again later.
  file.exceptions(std::ios::badbit | std::ios::failbit);

  // c_str explicitly converts to c String.
  String nname = name;

  if (nname.size() < 3 || nname.substr(nname.length() - 3, 3) != ".gz") {
    nname += ".gz";
  }

  try {
    file.open(nname.c_str());
  } catch (const std::ios::failure&) {
    std::ostringstream os;
    os << "Cannot open output file: " << nname << '\n'
       << "Maybe you don't have write access "
       << "to the directory or the file?";
    throw std::runtime_error(os.str());
  }

  // See if the file is ok.
  // FIXME: This should not be necessary anymore in the future, when
  // g++ stream exceptions work properly. (In that case we would not
  // get here if there really was a problem, because of the exception
  // thrown by open().)
  if (!file) {
    std::ostringstream os;
    os << "Cannot open output file: " << nname << '\n'
       << "Maybe you don't have write access "
       << "to the directory or the file?";
    throw std::runtime_error(os.str());
  }
}

#endif /* ENABLE_ZLIB */

//! Open file for XML input
/*!
  This function opens an XML file for reading.

  \param ifs   Input filestream
  \param name  Filename
*/
void xml_open_input_file(std::ifstream& ifs, const String& name) {
  // Tell the stream that it should throw exceptions.
  // Badbit means that the entire stream is corrupted.
  // On the other hand, end of file will not lead to an exception, you
  // have to check this manually!
  ifs.exceptions(std::ios::badbit);

  // c_str explicitly converts to c String.
  try {
    ifs.open(name.c_str());
  } catch (const std::ios::failure&) {
    std::ostringstream os;
    os << "Cannot open input file: " << name << '\n'
       << "Maybe the file does not exist?";
    throw std::runtime_error(os.str());
  }

  // See if the file is ok.
  // FIXME: This should not be necessary anymore in the future, when
  // g++ stream exceptions work properly.
  if (!ifs) {
    std::ostringstream os;
    os << "Cannot open input file: " << name << '\n'
       << "Maybe the file does not exist?";
    throw std::runtime_error(os.str());
  }
}

#ifdef ENABLE_ZLIB

//! Open file for zipped XML input
/*!
  This function opens a zipped XML file for reading.

  \param ifs   Input filestream
  \param name  Filename
*/
void xml_open_input_file(igzstream& ifs, const String& name) {
  // Tell the stream that it should throw exceptions.
  // Badbit means that the entire stream is corrupted.
  // On the other hand, end of file will not lead to an exception, you
  // have to check this manually!
  ifs.exceptions(std::ios::badbit);

  // c_str explicitly converts to c String.
  try {
    ifs.open(name.c_str());
  } catch (const std::ios::failure&) {
    std::ostringstream os;
    os << "Cannot open input file: " << name << '\n'
       << "Maybe the file does not exist?";
    throw std::runtime_error(os.str());
  }

  // See if the file is ok.
  // FIXME: This should not be necessary anymore in the future, when
  // g++ stream exceptions work properly.
  if (!ifs) {
    std::ostringstream os;
    os << "Cannot open input file: " << name << '\n'
       << "Maybe the file does not exist?";
    throw std::runtime_error(os.str());
  }
}

#endif /* ENABLE_ZLIB */

////////////////////////////////////////////////////////////////////////////
//   General XML functions (file header, start root tag, end root tag)
////////////////////////////////////////////////////////////////////////////

//! Throws XML parser runtime error
/*!
  This is used quite often inside the parsing routines so it's a
  function for itself.

  \param str_error Error description
*/
void xml_parse_error(const String& str_error) {
  std::ostringstream os;
  os << "XML parse error: " << str_error << '\n'
     << "Check syntax of XML file\n";
  throw std::runtime_error(os.str());
}

//! Throws XML parser runtime error
/*!
  This is used quite often inside the data parsing routines so it's a
  function for itself.

  \param tag        XMLTag
  \param str_error  Error description
*/
void xml_data_parse_error(XMLTag& tag, const String& str_error) {
  std::ostringstream os;
  os << "XML data parse error: Error reading ";
  tag.write_to_stream(os);
  os << str_error << "\n"
     << "Check syntax of XML file. A possible cause is that the file "
     << "contains NaN or Inf values.\n";
  throw std::runtime_error(os.str());
}

//! Reads XML header and root tag
/*!
  Check whether XML file has correct version tag and reads arts root
  tag information.

  \param is     Input stream
  \param ftype  File type
  \param ntype  Numeric type
  \param etype  Endian type
*/
void xml_read_header_from_stream(std::istream& is,
                                 FileType& ftype,
                                 NumericType& ntype,
                                 EndianType& etype) {
  char str[6];
  std::stringbuf strbuf;
  XMLTag tag;
  String strtype;

  while (!is.fail() && isspace(is.peek())) is.get();

  is.get(str, 6, ' ');

  if (std::string(str) != "<?xml") {
    xml_parse_error(
        "Input file is not a valid xml file "
        "(<?xml not found)");
  }

  is.get(strbuf, '>');
  is.get();

  if (is.fail()) {
    xml_parse_error("Input file is not a valid xml file");
  }

  tag.read_from_stream(is);
  tag.check_name("arts"sv);

  Index version{};
  tag.get_attribute_value("version", version);
  if (version > ARTS_XML_VERSION) {
    throw std::runtime_error(std::format(
        R"(Bad version in "arts"-tag.  Compiled version number: {}, got {})",
        ARTS_XML_VERSION,
        version));
  }

  // Check file format
  tag.get_attribute_value("format", strtype);
  if (strtype == "binary") {
    ftype = FileType::binary;
  } else {
    ftype = FileType::ascii;
  }

  // Check endian type
  tag.get_attribute_value("endian_type", strtype, "little");
  if (strtype == "little") {
    etype = ENDIAN_TYPE_LITTLE;
  } else if (strtype == "big") {
    etype = ENDIAN_TYPE_BIG;
  } else {
    std::ostringstream os;
    os << "  Error: Unknown endian type \"" << strtype
       << "\" specified in XML file.\n";
    throw std::runtime_error(os.str());
  }

  // Check numeric type
  tag.get_attribute_value("numeric_type", strtype, "double");
  if (strtype == "float") {
    ntype = NUMERIC_TYPE_FLOAT;
  } else if (strtype == "double") {
    ntype = NUMERIC_TYPE_DOUBLE;
  } else {
    std::ostringstream os;
    os << "  Error: Unknown numeric type \"" << strtype
       << "\" specified in XML file.\n";
    throw std::runtime_error(os.str());
  }
}

//! Reads closing root tag
/*!
  Checks whether XML file ends correctly with \</arts\>.

  \param is  Input stream
*/
void xml_read_footer_from_stream(std::istream& is) {
  XMLTag tag;

  tag.read_from_stream(is);
  tag.check_name("/arts"sv);
}

//! Writes XML header and root tag
/*!
  \param os     Output stream
  \param ftype  File type
*/
void xml_write_header_to_stream(std::ostream& os, FileType ftype) {
  XMLTag tag;

  std::println(os, "<?xml version=\"1.0\"?>");

  tag.set_name("arts");
  switch (ftype) {
    case FileType::ascii:
    case FileType::zascii: tag.add_attribute("format", "ascii"); break;
    case FileType::binary: tag.add_attribute("format", "binary"); break;
  }

  tag.add_attribute("version", ARTS_XML_VERSION);

  tag.write_to_stream(os);

  std::println(os);
}

//! Write closing root tag
/*!
  \param os Output stream
*/
void xml_write_footer_to_stream(std::ostream& os) {
  XMLTag tag;

  tag.set_name("/arts");
  tag.write_to_stream(os);

  std::println(os);
}

void xml_set_stream_precision(std::ostream& os) {
  os << std::setprecision(std::numeric_limits<Numeric>::max_digits10);
}

//! Get the content of an xml tag as a string
void parse_xml_tag_content_as_string(std::istream& is_xml, String& content) {
  char dummy;

  content = "";
  dummy   = (char)is_xml.peek();
  while (is_xml && dummy != '<') {
    is_xml.get(dummy);
    content += dummy;
    dummy    = (char)is_xml.peek();
  }

  if (!is_xml) throw std::runtime_error("Unexpected end of file.");
}
