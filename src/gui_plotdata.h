#ifndef PLOTDATA_H
#define PLOTDATA_H

#include "../3rdparty/gui/imgui/imgui.h"
#include <cmath>
#include <functional>
#include <memory>
#include <string>
#include <vector>

template<class T>
class PlotData {
  std::string mname;
  T one;
  T two;
  bool cur;
public:
  PlotData(const std::string& name, std::size_t n) : mname(name), one(n), two(n), cur(true) {}
  const T& get() const {if(cur) return one; else return two;}
  bool set(const T& x) {if (x.nelem() not_eq size()) return false; if(cur) two = x; else one = x; cur = not cur; return true;}
  const std::string& name() const {return mname;}
  std::size_t size() const {return std::size_t(one.nelem());}
};

enum class LineType {
  Average,
};

typedef ImVec2 (* LineGetter)(void *, int);

// Change in C++20 to use the atmopic locks
template <class T>
class Line {
  std::string mname;
  PlotData<T>* mx;
  PlotData<T>* my;
  bool xy;
  LineType mtype;
  size_t typemodifier;
  
  float x(int i) {if (xy) return mx -> get()[i]; else return i;}
  float y(int i) {return my -> get()[i];}
  
  float avg_x(int i0) {float sum_x=0.0f; for(size_t i=typemodifier*i0; i<typemodifier*i0+typemodifier; i++) sum_x += x(i); return sum_x / typemodifier;}
  float avg_y(int i0) {float sum_y=0.0f; for(size_t i=typemodifier*i0; i<typemodifier*i0+typemodifier; i++) sum_y += y(i); return sum_y / typemodifier;}
  LineGetter avg() const
  {
    return [](void * data, int i) {
      return ImVec2(static_cast<Line<T> *>(data) -> avg_x(i),
                    static_cast<Line<T> *>(data) -> avg_y(i));
    };
  }
  
public:
  Line(const std::string& name, PlotData<T>* x, PlotData<T>* y) : mname(name), mx(x), my(y), xy(true), mtype(LineType::Average), typemodifier(1) {}
  Line(const std::string& name, PlotData<T>* y) : mname(name), my(y), xy(false), mtype(LineType::Average), typemodifier(1) {}
  const char * name() const {return mname.c_str();}
  
  void changeX(PlotData<T>* x) {mx = x; xy=true;}
  void changeY(PlotData<T>* y) {my = y;}
  void swapaxis() {std::swap(mx, my);}
  
  // LineType data
  int size() const
  {
    switch (mtype) {
      case LineType::Average:
        return int(my -> size()) / typemodifier;
    }
    return 0;
  }
  
  LineGetter getter() const
  {
    switch (mtype) {
      case LineType::Average: return avg();
    }
    return [](void*, int){return ImVec2();};
  }
  
  // LineType setter
  size_t maxsize() const {return my -> size();}
  size_t TypeModifier() const {return typemodifier;}
  void Linear(size_t n) {mtype = LineType::Average; if(n) typemodifier=n; else typemodifier=1;}
};

#endif  // PLOTDATA_H
