// Downloads or reuses precise SP3 orbit files for several UTC epochs and prints
// the discovered satellite coverage. Requires Earthdata/CDDIS access to be
// configured as described in the docs.
#include <lupnt/lupnt.h>

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <vector>

using namespace lupnt;

int main() {
  const std::vector<Real> epochs_utc = {
      GregorianToTime(2025, 1, 1, 0, 0, 0.0),
      GregorianToTime(2025, 1, 2, 0, 0, 0.0),
      GregorianToTime(2025, 1, 3, 12, 0, 0.0),
  };

  for (const Real epoch_utc : epochs_utc) {
    std::cout << "\n=== " << TimeToGregorianString(epoch_utc, 0) << " UTC ===\n";
    std::cout << "Expected SP3: " << Sp3Loader::FilenameForEpoch(epoch_utc, Time::UTC) << "\n";
    std::cout << "CDDIS URL: " << Sp3Loader::UrlForEpoch(epoch_utc, Time::UTC) << "\n";

    const std::filesystem::path filepath = Sp3Loader::DownloadFileForEpoch(epoch_utc, Time::UTC);
    std::cout << "Downloaded/loaded file:\n  " << filepath.string() << "\n";

    Sp3Loader loader(filepath);
    const auto& sats = loader.GetSatellites();
    std::cout << "Number of satellites: " << sats.size() << "\n";
    std::cout << "First satellites:";
    for (size_t i = 0; i < std::min<size_t>(8, sats.size()); ++i) {
      std::cout << (i == 0 ? " " : ", ") << sats[i];
    }
    std::cout << "\n";

    if (!sats.empty()) {
      auto [t_min, t_max] = loader.GetTimeSpan(sats.front());
      std::cout << "Epoch coverage for " << sats.front() << ": " << t_min << " to " << t_max
                << " TAI seconds\n";
    }
  }

  return 0;
}
