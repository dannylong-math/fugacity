//
// Tests for synthesize::EoS construction preconditions.
//
// The EoS constructor requires the ideal and residual contributions to describe
// the same number of components. In a debug build that precondition is checked
// by SYNTHESIZE_ASSERT, which throws std::logic_error on a mismatch; in a release
// build the check is elided (so the throwing test is compiled only under debug).
//
#include "support/analytic_eos_models.hpp"
#include "synthesize/core/eos_pair.hpp"

#include <boost/ut.hpp>
#include <stdexcept>

using namespace boost::ut;
using namespace synthesize_test;

int main()
{
    suite<"eos_pair"> s = [] {
        // A well-formed pair (matching component counts) must construct cleanly.
        "matched component counts construct"_test = [] {
            auto eos = make_binary_model();
            expect(eos.size() == 2_ul);
        };

#ifndef NDEBUG
        // Mismatched component counts: the ideal model describes 2 components and
        // the residual describes 3, so the constructor precondition must reject
        // the pair by throwing std::logic_error (SYNTHESIZE_ASSERT, debug only).
        "mismatched component counts throw std::logic_error"_test = [] {
            IdealGasTestModel<2> ideal{{2.5, 3.1}, {1.5, 2.0}};
            VirialResidualTestModel<3> residual{
                {1.0e-3, 1.5e-3, 1.0e-3}, {1.0e-1, 8.0e-2, 5.0e-2}, {2.0e-6, 1.0e-6, 1.0e-6}};
            expect(throws<std::logic_error>([&] {
                const synthesize::EoS<IdealGasTestModel<2>, VirialResidualTestModel<3>> bad{ideal, residual};
                (void)bad;
            }));
        };
#endif
    };

    return ::boost::ut::cfg<>.run();
}
