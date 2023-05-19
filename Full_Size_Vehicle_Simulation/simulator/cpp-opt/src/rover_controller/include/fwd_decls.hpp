namespace roahm {
struct CostFcnInfo;
struct CudaInfo;
class DynObs;
struct FlZonoObsSet;
struct ZonoInfo;
struct Constraints;
class Gencon2Out;
struct FrsMega;
struct FrsSelectInfo;
struct IndividualMuSigma;
enum class ManuType;
struct PointXY;
struct PointXYH;
struct ProbIntegrationInputs;
struct RiskProblemDescription;
struct MuSigmaMulti;
struct MaybeRisky;
struct RiskRtdPlanningOutputs;
class RiskRtd;
class RoverState;
struct Sliceable;
struct SlicedInfo;
struct ZonoSliceInfo;
struct Vehrs;
struct WaypointCost;
struct WaypointCostRet;

namespace risk_rtd_ipopt_problem {
struct OptimizationInputs;
class RiskRtdIpoptProblem;
} // namespace risk_rtd_ipopt_problem

namespace risk_cost {
struct RiskCostValues;
struct RiskCost;
} // namespace risk_cost
namespace risk_constraint {
struct RiskConstraintValues;
struct RiskConstraint;
} // namespace risk_constraint

namespace pre_slice {
struct PreSliceOutputs;
}

namespace fl_zono_constraint {
class FlZonoConstraint;
}

namespace fl_zono_ipopt_problem {
class FlZonoIpoptProblem;
}
} // namespace roahm