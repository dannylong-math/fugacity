#include "support/analytic_eos_models.hpp"
#include "fugacity/core/concepts.hpp"

#include <boost/ut.hpp>

using namespace boost::ut;
using namespace fugacity_test;
namespace fug = fugacity;

// A type that does not model EquationOfState at all.
namespace {
struct NotAnEoS {
    static int foo() { return 0; }
};
} // namespace

// Compile-time verification of the concept hierarchy. If any of these fire, the
// concepts in concepts.hpp (or the test models) regressed.
static_assert(fug::EquationOfState<IdealGasTestModel<2>>);
static_assert(fug::EquationOfState<VirialResidualTestModel<2>>);

static_assert(fug::IdealEoS<IdealGasTestModel<2>>);
static_assert(!fug::IdealEoS<VirialResidualTestModel<2>>); // not derived from BaseIdealEoS

static_assert(fug::ResidualEoS<VirialResidualTestModel<2>>);
static_assert(!fug::ResidualEoS<IdealGasTestModel<2>>); // it IS an ideal EoS

static_assert(!fug::EquationOfState<NotAnEoS>);
static_assert(!fug::IdealEoS<NotAnEoS>);
static_assert(!fug::ResidualEoS<NotAnEoS>);

// The dynamic-extent variants also model the concepts.
static_assert(fug::IdealEoS<IdealGasTestModel<3>>);
static_assert(fug::ResidualEoS<VirialResidualTestModel<3>>);

int main()
{
    suite<"concepts"> s = [] {
        "models satisfy the expected concepts"_test = [] {
            // The static_asserts above already enforce this at compile time;
            // mirror them at runtime so the test reports as executed.
            expect(fug::IdealEoS<IdealGasTestModel<2>>);
            expect(fug::ResidualEoS<VirialResidualTestModel<2>>);
            expect(not fug::EquationOfState<NotAnEoS>);
        };

        "EoS pair exposes size and members"_test = [] {
            auto eos = make_binary_model();
            expect(eos.size() == 2_ul);
        };
    };

    return ::boost::ut::cfg<>.run();
}
