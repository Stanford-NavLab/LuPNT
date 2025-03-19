// pybind11
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace lupnt;

class Animal {
public:
  virtual ~Animal() = default;
  virtual std::string go(int n_times) = 0;
  virtual std::string name() { return "unknown"; }
};
class Dog : public Animal {
public:
  std::string go(int n_times) override {
    std::string result;
    for (int i = 0; i < n_times; ++i) result += bark() + " ";
    return result;
  }
  virtual std::string bark() { return "woof!"; }
};
class Husky : public Dog {};

template <class AnimalBase = Animal> class PyAnimal : public AnimalBase {
public:
  using AnimalBase::AnimalBase;  // Inherit constructors
  std::string go(int n_times) override {
    PYBIND11_OVERRIDE_PURE(std::string, AnimalBase, go, n_times);
  }
  std::string name() override { PYBIND11_OVERRIDE(std::string, AnimalBase, name, ); }
};

template <class DogBase = Dog> class PyDog : public PyAnimal<DogBase> {
public:
  using PyAnimal<DogBase>::PyAnimal;  // Inherit constructors
  // Override PyAnimal's pure virtual go() with a non-pure one:
  std::string go(int n_times) override { PYBIND11_OVERRIDE(std::string, DogBase, go, n_times); }
  std::string bark() override { PYBIND11_OVERRIDE(std::string, DogBase, bark, ); }
};

void init_sandbox(py::module &m) {
  py::class_<Animal, PyAnimal<>> animal(m, "Animal");
  py::class_<Dog, Animal, PyDog<>> dog(m, "Dog");
  py::class_<Husky, Dog, PyDog<Husky>> husky(m, "Husky");

  // Add animal methods
  animal.def(py::init<>()).def("go", &Animal::go).def("name", &Animal::name);
  dog.def(py::init<>()).def("bark", &Dog::bark);
  husky.def(py::init<>());
}
