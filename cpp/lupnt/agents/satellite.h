#include "lupnt/agents/agent.h"

namespace lupnt {

  // Satellite
  class Satellite : public AgentWithDynamics {
  public:
    Satellite() = default;
    Satellite(Config& config);
    virtual void Log(Real time) override;
    virtual void LogCesium() override;
  };

}  // namespace lupnt
