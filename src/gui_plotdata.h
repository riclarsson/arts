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

class PlotData {
  String mname;
  Vector one;
  Vector two;
  bool cur;
public:
  PlotData(const String& name, std::size_t n) : mname(name), one(n), two(n), cur(true) {}
  const Vector& get() const {if(cur) return one; else return two;}
  bool set(const Vector& x) {if (x.nelem() not_eq nelem()) return false; if(cur) two = x; else one = x; cur = not cur; return true;}
  const String& name() const {return mname;}
  Index nelem() const {return one.nelem();}
};

enum class LineType {
  Average,
  RunningAverage
};

typedef ImVec2 (* LineGetter)(void *, int);

// Change in C++20 to use the atmopic locks
class Line {
  String mname;
  PlotData* mx;
  PlotData* my;
  bool xy;
  LineType mtype;
  std::size_t typemodifier;
  
  float x(int i) {if (xy) return float(mx -> get()[i]); else return float(i);}
  float y(int i) {return float(my -> get()[i]);}
  
  float avg_x(int i0) {float sum_x=0.0f; for(std::size_t i=typemodifier*i0; i<typemodifier*i0+typemodifier; i++) sum_x += x(int(i)); return sum_x / float(typemodifier);}
  float avg_y(int i0) {float sum_y=0.0f; for(std::size_t i=typemodifier*i0; i<typemodifier*i0+typemodifier; i++) sum_y += y(int(i)); return sum_y / float(typemodifier);}
  LineGetter avg() const
  {
    return [](void * data, int i) {
      return ImVec2(static_cast<Line *>(data) -> avg_x(i),
                    static_cast<Line *>(data) -> avg_y(i));
    };
  }
  
  float runavg_x(int i0) {float sum_x=0.0f; for(std::size_t i=i0; i<i0+typemodifier; i++) sum_x += x(int(i)); return sum_x / float(typemodifier);}
  float runavg_y(int i0) {float sum_y=0.0f; for(std::size_t i=i0; i<i0+typemodifier; i++) sum_y += y(int(i)); return sum_y / float(typemodifier);}
  LineGetter runavg() const
  {
    return [](void * data, int i) {
      return ImVec2(static_cast<Line *>(data) -> runavg_x(i),
                    static_cast<Line *>(data) -> runavg_y(i));
    };
  }
  
public:
  Line(const String& name, PlotData* x, PlotData* y) : mname(name), mx(x), my(y), xy(true), mtype(LineType::Average), typemodifier(1) {}
  Line(const String& name, PlotData* y) : mname(name), my(y), xy(false), mtype(LineType::Average), typemodifier(1) {}
  const char * name() const {return mname.c_str();}
  
  void changeX(PlotData* x) {mx = x; xy=true;}
  void changeY(PlotData* y) {my = y;}
  void swapaxis() {std::swap(mx, my);}
  
  // LineType data
  int size() const
  {
    switch (mtype) {
      case LineType::Average:
        return int(my -> nelem() / typemodifier);
      case LineType::RunningAverage:
        return int(my -> nelem() - typemodifier) + 1;
    }
    return 0;
  }
  
  LineGetter getter() const
  {
    switch (mtype) {
      case LineType::Average: return avg();
      case LineType::RunningAverage: return runavg();
    }
    return [](void*, int){return ImVec2();};
  }
  
  std::size_t maxsize() const
  {
    switch (mtype) {
      case LineType::Average:
        return my -> nelem() / 2;
      case LineType::RunningAverage:
        return my -> nelem() / 2;
    }
    return 0;
  }
  
  Index TypeModifier() const {return typemodifier;}
  void Linear(std::size_t n) {mtype = LineType::Average; if(n) typemodifier=n; else if(n > maxsize()) typemodifier = maxsize(); else typemodifier=1;}
  void RunningLinear(std::size_t n) {mtype = LineType::RunningAverage; if(n) typemodifier=n; else if(n > maxsize()) typemodifier = maxsize(); else typemodifier=1;}
};

#endif  // GUI_PLOTDATA_H
