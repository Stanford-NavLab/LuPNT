#include <lupnt/core/scheduler.h>

#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <queue>
#include <vector>

using namespace std;
using namespace lupnt;

// Transmission class
struct Transmission {
  int id_;
};

class Transceiver;

// Channel class
class Channel {
private:
  int id_;
  double delay_ = 0.5;
  vector<shared_ptr<Transceiver>> transceivers_;

public:
  void Send(double t, const Transmission &transm, const Transceiver &sender);
  void Add(shared_ptr<Transceiver> transceiver) { transceivers_.push_back(transceiver); }
  vector<shared_ptr<Transceiver>> GetTransceivers() { return transceivers_; }
};

// Transceiver class
class Transceiver {
private:
  shared_ptr<Channel> channel_;
  function<void(double, const Transmission &)> receive_callback_;

public:
  Transceiver(shared_ptr<Channel> channel) : channel_(channel) {}
  void SetReceiveCallback(function<void(double, const Transmission &)> callback) {
    receive_callback_ = callback;
  }

  void Send(double t, const Transmission &transm) { channel_->Send(t, transm, *this); }

  void Receive(double t, const Transmission &transm) { receive_callback_(t, transm); }
};

void Channel::Send(double t, const Transmission &transm, const Transceiver &sender) {
  for (const auto &transceiver : transceivers_) {
    if (transceiver.get() != &sender) {
      Scheduler::Schedule(t + delay_,
                          [transceiver, transm](double t) { transceiver->Receive(t, transm); });
    }
  }
}

// Agent class
class RingAgent {
private:
  static int id_counter_;
  const int id_;
  vector<shared_ptr<Transceiver>> transceivers_;

public:
  RingAgent() : id_(id_counter_++), transceivers_(vector<shared_ptr<Transceiver>>()) {}
  int GetId() { return id_; }

  void Add(shared_ptr<Transceiver> transceiver) { transceivers_.push_back(transceiver); }
  vector<shared_ptr<Transceiver>> GetTransceivers() { return transceivers_; }
};

int RingAgent::id_counter_ = 0;

// Application class
class LeaderElectionSyncRingApp : public Application {
private:
  shared_ptr<RingAgent> agent_;
  int id_;
  int id_received_ = -1;
  bool is_leader_ = false;

public:
  LeaderElectionSyncRingApp(shared_ptr<RingAgent> agent) : agent_(agent), id_(agent->GetId()) {}
  double GetFrequency() override { return 1.0; }
  void Setup() override {}
  void Step(double t) override {
    if (id_received_ == -1) {
      cout << "[Agent " << id_ << "] Sending " << id_ << " at t = " << t << endl;
      for (const auto &transceiver : agent_->GetTransceivers()) {
        // transmittion instance with _id
        transceiver->Send(t, Transmission{id_});
      }
    } else if (id_received_ > id_) {
      cout << "[Agent " << id_ << "] Sending " << id_received_ << " at t = " << t << endl;
      for (const auto &transceiver : agent_->GetTransceivers()) {
        transceiver->Send(t, Transmission{id_received_});
      }
    } else if (id_received_ < id_) {
      cout << "[Agent " << id_ << "] Doing nothing" << " at t = " << t << endl;
    } else {
      is_leader_ = true;
      cout << "[Agent " << id_ << "] Leader" << " at t = " << t << endl;
    }
  }
  void TransmissionReceived(double t, const Transmission &transm) {
    cout << "[Agent " << id_ << "] Received " << transm.id_ << " at t = " << t << endl;
    id_received_ = max(id_received_, transm.id_);
  }
};

int main() {
  int n = 3;
  vector<shared_ptr<RingAgent>> agents(n);
  vector<shared_ptr<LeaderElectionSyncRingApp>> apps(n);
  vector<shared_ptr<Channel>> channels(n);

  // Create agents, applications, and channels
  for (int i = 0; i < n; i++) {
    agents[i] = make_shared<RingAgent>();
    apps[i] = make_shared<LeaderElectionSyncRingApp>(agents[i]);
    channels[i] = make_shared<Channel>();
  }

  // Create transceivers and add them to agents and channels
  for (int i = 0; i < n; i++) {
    auto ch1 = channels[i];
    auto ch2 = channels[(i + 1) % n];
    auto transc1 = make_shared<Transceiver>(ch1);
    auto transc2 = make_shared<Transceiver>(ch2);
    transc1->SetReceiveCallback(bind(&LeaderElectionSyncRingApp::TransmissionReceived,
                                     apps[i].get(), placeholders::_1, placeholders::_2));
    transc2->SetReceiveCallback(bind(&LeaderElectionSyncRingApp::TransmissionReceived,
                                     apps[i].get(), placeholders::_1, placeholders::_2));
    agents[i]->Add(transc1);
    agents[i]->Add(transc2);
    ch1->Add(transc1);
    ch2->Add(transc2);
  }

  // Schedule applications
  double t_start = 0.0;
  double freq = 1.0;
  for (int i = 0; i < n; i++) {
    Scheduler::ScheduleApplication(*apps[i], t_start, freq);
  }

  // Run simulation
  Scheduler::RunSimulation(4.0);
  return 0;
}
