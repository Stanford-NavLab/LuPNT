#include "lupnt/agents/agent.h"

namespace lupnt {

  // Rover
  class Rover : public AgentWithDynamics {
  protected:
    Cart3 r_ref_pa_;
    Cart6 rv_enu_;
    SurfaceState2D state2d_;

  public:
    Rover(Config& agent_config);
    virtual void Log(Real time) override;
    virtual void Propagate(Real t) override;
  };

}  // namespace lupnt
