#pragma once

#include <unordered_set>

#include "lupnt/core/config.h"
#include "lupnt/core/definitions.h"

namespace lupnt {

  class Device;

  class Channel {
  public:
    Channel() = default;

    /// @brief Construct a channel from a YAML config node, reading its `name`
    /// field (defaults to "channel" if absent).
    ///
    /// Invoked by Channel::FromConfig / AssetFactory when building the
    /// simulation's communication topology from a scenario config file.
    Channel(Config& config);
    ~Channel() = default;

    /// @brief Register a device as a participant on this channel.
    ///
    /// Called when wiring up agents/devices during simulation setup so that
    /// `Send`/`Receive` know which other devices share this channel.
    void AddDevice(Device* device);

    /// @brief Name of this channel (as set from config, or "channel" by default).
    std::string GetName() const { return name_; }

    /// @brief Broadcast `data` from transmitting device `tx` to every other
    /// registered device on the channel that can receive (a `Receiver` or
    /// `Transponder`), invoking their `Receive(t, data)`.
    ///
    /// Called by a `Transmitter`/`Transponder` device when it sends data at
    /// time `t`, to deliver that data to all other devices sharing this
    /// channel (e.g. a ground-station receiver picking up a satellite's
    /// transmission).
    ///
    /// @param tx   Sending device (must be non-null; skipped as its own recipient)
    /// @param t    Simulation time of the send [s]
    /// @param data Opaque pointer to the data payload being sent (must be non-null)
    void Send(Device* tx, Real t, void* data);

    /// @brief Trigger every other registered device on the channel that can
    /// transmit (a `Transmitter` or `Transponder`) to send at time `t`,
    /// causing receiving device `rx` to get their data via `Send`.
    ///
    /// Called by a `Receiver`/`Transponder` device when it samples/listens at
    /// time `t`, to pull data from all transmitting devices sharing this channel.
    ///
    /// @param rx Receiving device (must be non-null; skipped as its own sender)
    /// @param t  Simulation time of the receive [s]
    void Receive(Device* rx, Real t);

    /// @brief Construct a `Channel` (or registered subclass) from a YAML
    /// config node's `class` field, via the `AssetFactory<Channel, Config&>` registry.
    static std::shared_ptr<Channel> FromConfig(Config& config);

  protected:
    std::string name_;
    std::unordered_set<Device*> devices_;
  };

}  // namespace lupnt
