// Author: Hongliang Gao;   Created: Jan 27 2011
// Description: sample predictors for cbp3.

#include <stdio.h>
#include <cassert>
#include <string.h>
#include <inttypes.h>

using namespace std;
#include "cbp3_def.h"
#include "cbp3_framework.h"


void PredictorInit() {
}

void PredictorReset() {
    // this function is called before EVERY run
    // it is used to reset predictors and change configurations

}

void PredictorRunACycle() {
    // get info about what uops are processed at each pipeline stage
    const cbp3_cycle_activity_t *cycle_info = get_cycle_info();

    // make prediction at fetch stage
    for (int i = 0; i < cycle_info->num_fetch; i++) {
        uint32_t fe_ptr = cycle_info->fetch_q[i];
        const cbp3_uop_dynamic_t *uop = &fetch_entry(fe_ptr)->uop;

    }

    for (int i = 0; i < cycle_info->num_retire; i++) {
        uint32_t rob_ptr = cycle_info->retire_q[i];
        const cbp3_uop_dynamic_t *uop = &rob_entry(rob_ptr)->uop;

    }
}

void PredictorRunEnd() {
}

void PredictorExit() {
}
