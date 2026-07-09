#pragma once

#include "lupnt/core/definitions.h"
#include "lupnt/devices/device.h"
#include "lupnt/measurements/channel.h"

namespace lupnt {

  /// @brief Generic radio transmitter device: binds to a named `Channel` and
  /// pushes outgoing data to it.
  ///
  /// Used as the base for link-budget/measurement-generating transmitters
  /// (e.g. `GnssTransmitter` in `space_comms.h`) attached to a satellite or
  /// ground-station `Agent`; the bound `Channel` (resolved in `Setup`) is the
  /// shared medium used by `measurements/` models to compute link
  /// observables between a transmitter and a `Receiver`/`Transponder`.
  class Transmitter : public Device {
  public:
    Transmitter() = default;
    /// @brief Construct a transmitter from a YAML config node (see
    /// `Device::Device(Config&)`).
    Transmitter(Config& config);
    /// @brief Transmitter device step; currently a no-op placeholder.
    /// Overrides `Device::Step`.
    virtual void Step(Real t) override;
    /// @brief Resolve and bind `channel_` from `config_["channel"]` via the
    /// owning agent's `Simulation::GetChannel`. Overrides `Device::Setup`;
    /// logs a warning if no `"channel"` key is configured.
    virtual void Setup() override;

    /// @brief Notify the bound `channel_` that this transmitter is ready to
    /// send, without an explicit payload (the channel pulls data via
    /// `AddData`/`data_`). Requires `channel_` to be set (see `Setup`).
    ///
    /// @param t  Simulation time of the send event [s]
    virtual void Send(Real t);

    /// @brief Push `data` onto the bound `channel_` at time `t`. Requires
    /// `channel_` and `data` to be non-null (see `Setup`).
    ///
    /// @param t     Simulation time of the send event [s]
    /// @param data  Pointer to the payload to send (ownership/type defined by
    ///              the specific `Channel`/measurement model)
    virtual void Send(Real t, void* data);

    /// @brief Queue `data` for the next `Send` call.
    virtual void AddData(void* data) { data_.push_back(data); }
    /// @brief Clear all queued outgoing data.
    virtual void EmptyData() { data_.clear(); }

  protected:
    Channel* channel_;
    std::vector<void*> data_;
  };

  /// @brief Generic radio receiver device: binds to a named `Channel` and
  /// collects incoming data from it.
  ///
  /// Used as the base for measurement-generating receivers (e.g.
  /// `GnssReceiver` in `space_comms.h`) attached to a satellite/rover/
  /// ground-station `Agent`; received payloads accumulate in `data_` for
  /// downstream `measurements/` processing.
  class Receiver : public Device {
  public:
    Receiver() = default;
    /// @brief Construct a receiver from a YAML config node (see
    /// `Device::Device(Config&)`).
    Receiver(Config& config);
    /// @brief Receiver device step; currently a no-op placeholder. Overrides
    /// `Device::Step`.
    virtual void Step(Real t) override;
    /// @brief Resolve and bind `channel_` from `config_["channel"]` via the
    /// owning agent's `Simulation::GetChannel`. Overrides `Device::Setup`;
    /// logs a warning if no `"channel"` key is configured.
    virtual void Setup() override;

    /// @brief Notify this receiver that the channel has data available,
    /// without an explicit payload.
    ///
    /// @param t  Simulation time of the receive event [s]
    virtual void Receive(Real t);

    /// @brief Append an incoming payload `data` (received via `channel_`) to
    /// `data_` at time `t`.
    ///
    /// @param t     Simulation time of the receive event [s]
    /// @param data  Pointer to the received payload (ownership/type defined
    ///              by the specific `Channel`/measurement model)
    virtual void Receive(Real t, void* data);

  protected:
    Channel* channel_;
    std::vector<void*> data_;
  };

  /// @brief Combined transmit/receive radio device: binds to a named
  /// `Channel` and both sends and receives data through it.
  ///
  /// Models a transponder-type radio link (e.g. for two-way range/Doppler
  /// measurements between agents) attached to a satellite/rover/
  /// ground-station `Agent`, combining the `Transmitter`/`Receiver` send and
  /// receive behaviors on a single `Channel`.
  class Transponder : public Device {
  public:
    Transponder() = default;
    /// @brief Construct a transponder from a YAML config node (see
    /// `Device::Device(Config&)`).
    Transponder(Config& config);
    /// @brief Transponder device step; currently a no-op placeholder.
    /// Overrides `Device::Step`.
    virtual void Step(Real t) override;
    /// @brief Resolve and bind `channel_` from `config_["channel"]` via the
    /// owning agent's `Simulation::GetChannel`. Overrides `Device::Setup`;
    /// logs a warning if no `"channel"` key is configured.
    virtual void Setup() override;

    /// @brief Notify the bound `channel_` that this transponder is ready to
    /// send, without an explicit payload. Requires `channel_` to be set (see
    /// `Setup`). Same role as `Transmitter::Send(Real)`.
    virtual void Send(Real t);

    /// @brief Push `data` onto the bound `channel_` at time `t`. Requires
    /// `channel_` and `data` to be non-null. Same role as
    /// `Transmitter::Send(Real, void*)`.
    virtual void Send(Real t, void* data);

    /// @brief Request data from the bound `channel_` at time `t` (calls
    /// `Channel::Receive`). Requires `channel_` to be set.
    ///
    /// @param t  Simulation time of the receive event [s]
    virtual void Receive(Real t);

    /// @brief Append an incoming payload `data` (provided by `channel_`) to
    /// `data_` at time `t`. Requires `data` to be non-null. Same role as
    /// `Receiver::Receive(Real, void*)`.
    virtual void Receive(Real t, void* data);

  protected:
    Channel* channel_;
    std::vector<void*> data_;
  };
}  // namespace lupnt
