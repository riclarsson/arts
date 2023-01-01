#include <exception>
#include <iostream>

#include "atm.h"
#include "gridded_fields.h"
#include "species.h"
#include "species_tags.h"

void point() {
  using namespace Atm;
  using enum Key;
  Point atm{pressure, 3e4, temperature, 250, ArrayOfSpeciesTag{"O2-66"}, 0.21};
  std::cout << atm << '\n' << '\n';

  atm.set(ArrayOfSpeciesTag{"H2O-161"}, 0.01);
  std::cout << atm << '\n' << '\n';

  std::cout << atm[ArrayOfSpeciesTag{"O2-66"}] << '\n' << '\n';
  std::cout << atm[pressure] << ' ' << atm.P() << '\n' << '\n';

  ARTS_USER_ERROR_IF( not atm.has(pressure), "should have pressure")
  ARTS_USER_ERROR_IF( not atm.has_data(wind_u, mag_v), "should have fields")
  ARTS_USER_ERROR_IF( not atm.has(ArrayOfSpeciesTag{"O2-66"}), "should have O2-66")
  ARTS_USER_ERROR_IF( atm.has(ArrayOfSpeciesTag{"N2-44"}), "should not have N2-44")

  std::cout << "Mean mass (only O2 and Water): " << atm.mean_mass() << '\n';
  atm.set(ArrayOfSpeciesTag{"N2-44"}, 0.78);
  std::cout << "Mean mass (added N2): " << atm.mean_mass() << '\n';
}

void field() {
  using namespace Atm;
  using enum Key;
  Field atm{temperature, 250.0, ArrayOfSpeciesTag{"O2-66"}, 0.21};
  atm.top_of_atmosphere(100e3);
  std::cout << atm << '\n' << '\n';

  atm.set(ArrayOfSpeciesTag{"H2O-161"}, 0.01);
  std::cout << atm << '\n' << '\n';

  GriddedField3 data;
  data.set_grid(0, {0, 2, 4});
  data.set_grid_name(0, "Pressure");
  data.set_grid(1, {1e-3, 2, 4});
  data.set_grid_name(1, "Latitude");
  data.set_grid(2, {0, 2, 4});
  data.set_grid_name(2, "Longitude");
  data.data.resize(3, 3, 3);
  for (Index i1 = 0; i1 < data.get_grid_size(0); i1++) {
    for (Index i2 = 0; i2 < data.get_grid_size(1); i2++) {
      for (Index i3 = 0; i3 < data.get_grid_size(2); i3++) {
        data.data(i1, i2, i3) = static_cast<Numeric>(1 + i1 + 2 * i2 + 3 * i3);
      }
    }
  }
  atm.set(mag_u, data);

  std::cout << atm.at(1e-2, 3.0, 0.0) << '\n' << '\n';

  std::cout << atm.regularize({1e-1, 1.0}, {0.5}, {0.5}) << '\n'
            << '\n';
}

int main() try {
  std::cout << "\n"
               "##################################################"
               "\n"
               "START OF POINT"
               "\n"
               "\n";
  point();
  std::cout << "\n"
               "END   OF POINT"
               "\n"
               "##################################################"
               "\n"
               "\n";

  std::cout << "\n"
               "##################################################"
               "\n"
               "START OF FIELD"
               "\n"
               "\n";
  field();
  std::cout << "\n"
               "END   OF FIELD"
               "\n"
               "##################################################"
               "\n"
               "\n";
} catch (std::exception& e) {
  std::cerr << e.what() << '\n';
}