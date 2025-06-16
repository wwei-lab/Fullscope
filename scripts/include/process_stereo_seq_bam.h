#ifndef HANDLESTEREOBAM_H
#define HANDLESTEREOBAM_H

#include "utils.h"

using namespace std;
namespace HandleStereoBam {
    void build_reftable(string gtf_file, string bam_file, string refcid, string outfile, int num_threads, string tag);
}

#endif