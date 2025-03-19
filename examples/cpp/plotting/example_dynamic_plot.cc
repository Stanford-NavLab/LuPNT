#include <matplot/matplot.h>

#include <chrono>
#include <thread>

int main() {
  using namespace matplot;

  // Initial data
  std::vector<double> x(10), y(10);
  for (int i = 0; i < 10; ++i) {
    x[i] = i;
    y[i] = i * i;
  }

  // Create a figure and plot
  auto fig = figure(true);
  auto ax = fig->add_axes();
  auto l = ax->plot(x, y);
  ax->xlabel("X-axis");
  ax->ylabel("Y-axis");
  ax->title("Updating Plot Example");

  // Show initial plot
  // fig->show();

  // Update data in a loop
  for (int i = 10; i < 20; ++i) {
    std::this_thread::sleep_for(std::chrono::seconds(1));

    // Update data
    x.push_back(i);
    y.push_back(i * i);

    // Update plot data
    l->x_data(x);
    l->y_data(y);

    // Redraw the figure
    fig->draw();
  }

  return 0;
}
