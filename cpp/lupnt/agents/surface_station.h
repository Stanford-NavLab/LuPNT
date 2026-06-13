#include "lupnt/agents/agent.h"

namespace lupnt {

  class SurfaceStation : public AgentWithDynamics {
  public:
    SurfaceStation(Config& agent_config);
    virtual void Log(Real time) override;
    virtual void LogCesium() override;
  };

}  // namespace lupnt
