#ifndef GUI_PLOTDATA_H
#define GUI_PLOTDATA_H

#include "gui_inc.h"

#include <cmath>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "matpackI.h"
#include "mystring.h"

namespace ARTSGUI {
namespace Plotting {

/** Holds the plotting data --- might be modified to allow updating data more easily */
class Data {
  Vector vec;
  Numeric s, o;
public:
  Data(std::size_t n=0) : vec(n), s(1.0), o(0.0) {}
  Data(const Vector& v) : vec(v), s(1.0), o(0.0) {}
  Data(Vector&& v) : vec(std::move(v)), s(1.0), o(0.0) {}
  Numeric get(int i) const {return vec[i] / s - o / s;}
  bool set(const Vector& x) {if (x.nelem() not_eq nelem()) return false; else vec = x; return true;}
  bool set(ConstVectorView x) {if (x.nelem() not_eq nelem()) return false; else vec = x; return true;}
  void overwrite(Vector&& x) {vec = std::move(x);}
  Index nelem() const {return vec.nelem();}
  VectorView view() {return vec;}
  void scale(Numeric x) {s=x;}
  void offset(Numeric x) {o=x;}
  Numeric scale() const {return s;}
  Numeric offset() const {return o;}
};  // Data

using ArrayOfData = Array<Data>;

enum class LineType {
  Average,
  MatrixMultiplier
};  // LineType

/** LineGetter is a function pointer required by the plotting tool */
typedef ImVec2 (* LineGetter)(void *, int);

/** Holds the data of a single line */
class Line {
  String mname;
  Data* mx;
  Data* my;
  LineType mtype;
  std::size_t typemodifier;
  Matrix mmul;
  
  float x(int i) {if (mx not_eq nullptr) return float(mx -> get(i)); else return float(i);}
  float y(int i) {if (my not_eq nullptr) return float(my -> get(i)); else return 0.0f;}
  
  float avg_x(int i0) {float sum_x=0.0f; for(std::size_t i=typemodifier*i0; i<typemodifier*i0+typemodifier; i++) sum_x += x(int(i)); return sum_x / float(typemodifier);}
  float avg_y(int i0) {float sum_y=0.0f; for(std::size_t i=typemodifier*i0; i<typemodifier*i0+typemodifier; i++) sum_y += y(int(i)); return sum_y / float(typemodifier);}
  LineGetter avg() const
  {
    return [](void * data, int i) {
      return ImVec2(static_cast<Line *>(data) -> avg_x(i),
                    static_cast<Line *>(data) -> avg_y(i));
    };
  }
  
  float mul_x(int i0) {if (mx not_eq nullptr) return float(mmul(i0, joker) * mx -> view()); else return float(i0); }
  float mul_y(int i0) {if (my not_eq nullptr) return float(mmul(i0, joker) * my -> view()); else return float(i0); }
  LineGetter mul() const
  {
    return [](void * data, int i) {
      return ImVec2(static_cast<Line *>(data) -> mul_x(i),
                    static_cast<Line *>(data) -> mul_y(i));
    };
  }
  
public:
  Line(const String& name, Data* x, Data* y) : mname(name), mx(x), my(y), mtype(LineType::Average), typemodifier(1), mmul(0, 0) {}
  Line(const String& name, Data* y) : mname(name), mx(nullptr), my(y), mtype(LineType::Average), typemodifier(1), mmul(0, 0) {}
  Line(const String& name="") : mname(name), mx(nullptr), my(nullptr), mtype(LineType::Average), typemodifier(1), mmul(0, 0) {}
  
  void changeX(Data* x) {mx = x;}
  void changeY(Data* y) {my = y;}
  void swapaxis() {std::swap(mx, my);}
  
  // LineType data
  int size() const
  {
    if(my not_eq nullptr) {
      switch (mtype) {
        case LineType::Average: return int(my -> nelem() / typemodifier);
        case LineType::MatrixMultiplier: return int(mmul.nrows());
      }
    }
    return 0;
  }
  
  LineGetter getter() const
  {
    if(my not_eq nullptr) {
      switch (mtype) {
        case LineType::Average: return avg();
        case LineType::MatrixMultiplier: return mul();
      }
    }
    return [](void*, int){return ImVec2({NAN, NAN});};
  }
  
  std::size_t maxsize() const
  {
    if(my not_eq nullptr)
      return my -> nelem();
    else
      return 0;
  }
  
  // Types of averages
  Index TypeModifier() const {return typemodifier;}
  void Linear(std::size_t n) {
    mtype = LineType::Average;
    if (n > maxsize())
      typemodifier = maxsize();
    else if (n)
      typemodifier=n;
    else
      typemodifier=1;
  }
  void RunningLinear(std::size_t n) {
    mtype = LineType::MatrixMultiplier;
    if (n == 0) {
      mmul.resize(maxsize(), maxsize());
      n = maxsize();
    } else if(n >= maxsize()) {
      mmul.resize(1, maxsize());
      n = maxsize();
    } else if(n) {
      mmul.resize(maxsize() - n + 1, maxsize());
    } else {
      mmul.resize(0, 0);
      n = 1;
    }
    
    typemodifier = n;
    std::size_t j=0;
    const Numeric x = 1.0 / Numeric(n);
    for (Index i=0; i<mmul.nrows(); i++) {
      for (std::size_t k=j; k < j+n and Index(k)<mmul.ncols(); k++) {
        mmul(i, k) = x;
      }
      j += 1;
    }
  }
  
  // Names
  const String& name() const {return mname;}
  void name(const String& name) {mname = name;}
  
  // Dataviews
  Data * x() noexcept {return mx;}
  Data * y() noexcept {return my;}
};  // Line

using ArrayOfLine = Array<Line>;

class Frame {
  String mtitle;
  String mxlabel;
  String mylabel;
  ImPlotLimits mlimits;
  ArrayOfLine mlines;
public:
  Frame(const String& title, const String& x, const String& y) : mtitle(title), mxlabel(x), mylabel(y), mlimits(), mlines(0) {}
  Frame(const String& title, const String& x, const String& y, const Line& line) : mtitle(title), mxlabel(x), mylabel(y), mlimits(), mlines(1, line) {}
  Frame(const String& title, const String& x, const String& y, const ArrayOfLine& lines) : mtitle(title), mxlabel(x), mylabel(y), mlimits(), mlines(lines) {}
  
  // Names
  const String& title() const {return mtitle;}
  const String& xlabel() const {return mxlabel;}
  const String& ylabel() const {return mylabel;}
  void title(const String& x) {mtitle=x;}
  void xlabel(const String& x) {mxlabel=x;}
  void ylabel(const String& x) {mylabel=x;}
  
  // Line manipulation
  void push_back(Line& line) {mlines.push_back(line);}
  void pop_back() {mlines.pop_back();}
  void set_empty() {mlines.resize(0);}
  Line& operator[](Index i) {return mlines[i];}
  const Line& operator[](Index i) const {return mlines[i];}
  void lines(const ArrayOfLine& x) {mlines=x;}
  Index nelem() const {return mlines.nelem();}
  decltype(mlines.begin()) begin() {return mlines.begin();}
  decltype(mlines.end()) end() {return mlines.end();}
  decltype(mlines.cbegin()) cbegin() const {return mlines.cbegin();}
  decltype(mlines.cend()) cend() const {return mlines.cend();}
  
  // Range of frame
  ImPlotLimits limits() const {return mlimits;}
  void limits(ImPlotLimits x) {mlimits=x;}
  bool invalid_range() const {using std::isnan; return isnan(mlimits.X.Max) or isnan(mlimits.X.Min) or isnan(mlimits.Y.Max) or isnan(mlimits.Y.Min);}
};  // Frame

/** Plots a frame  in the current window */
bool PlotFrame(Frame& frame, Config& cfg, bool x_is_frequency);

};  // Plotting
};  // ARTSGUI

#endif  // GUI_PLOTDATA_H
