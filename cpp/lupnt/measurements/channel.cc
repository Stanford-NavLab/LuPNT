#include "lupnt/measurements/channel.h"

#include <string>

#include "lupnt/core/asset_factory.h"
#include "lupnt/core/logger.h"
#include "lupnt/devices/comms.h"
#include "lupnt/devices/device.h"

namespace lupnt {

  Channel::Channel(Config& config) {
    name_ = config["name"] ? config["name"].as<std::string>() : "channel";
    Logger::Debug(fmt::format("Creating {}", name_), "Channel");
  }

  void Channel::AddDevice(Device* device) { devices_.insert(device); }

  // Called by device sending data
  void Channel::Send(Device* tx, Real t, void* data) {
    LUPNT_CHECK(tx, "Device is null", "Channel");
    LUPNT_CHECK(data, "Data is null", "Channel");

    // Check if the sending device is a Transmitter or Transponder
    for (auto& rx : devices_) {
      if (rx == tx) continue;

      auto receiver = dynamic_cast<Receiver*>(rx);
      auto transponder = dynamic_cast<Transponder*>(rx);
      if (!(receiver || transponder)) continue;

      // Call Send on all other devices that can receive
      if (receiver) {
        receiver->Receive(t, data);
      } else if (transponder) {
        transponder->Receive(t, data);
      }
    }
  }

  // Called by device receiving data
  void Channel::Receive(Device* rx, Real t) {
    LUPNT_CHECK(rx, "Device is null", "Channel");

    for (auto& tx : devices_) {
      if (tx == rx) continue;  // Don't receive from self

      auto transmitter = dynamic_cast<Transmitter*>(tx);
      auto transponder = dynamic_cast<Transponder*>(tx);
      if (!(transmitter || transponder)) continue;

      // Call Send on all other devices that can send data
      if (transmitter) {
        transmitter->Send(t);
      } else if (transponder) {
        transponder->Send(t);
      }
    }
  }

  std::shared_ptr<Channel> Channel::FromConfig(Config& config) {
    std::string channel_class = config["class"].as<std::string>();
    Logger::Debug(fmt::format("Creating channel {}", channel_class), "Channel");
    return AssetFactory<Channel, Config&>::Create(channel_class, config);
  }

  REGISTER_FACTORY_CLASS(Channel, Channel);

  // Define the GetRegistry function for this specialization (must come before explicit
  // instantiation)
  template <> std::unordered_map<std::string, AssetFactory<Channel, Config&>::Creator>&
  AssetFactory<Channel, Config&>::GetRegistry() {
    return Registry();
  }

  // Explicit template instantiation to ensure single registry across library boundaries
  template class AssetFactory<Channel, Config&>;

}  // namespace lupnt
