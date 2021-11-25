#pragma once

#include "VectorField.hpp"
#include "Constants.hpp"

#include <string>
#include <fstream>

using namespace std;

void writeFile(string name, VectorField *m, float t)
{
    ofstream file;
    file.open(name);

    // Write out the header
    file << "# OOMMF OVF 2.0" << endl;
    file << "# Segment count: 1" << endl;
    file << "# Begin: Segment" << endl;
    file << "# Begin: Header" << endl;
    file << "# Title: m" << endl;
    file << "# meshtype: rectangular" << endl;
    file << "# meshunit: dimensionless" << endl;
    file << "# valueunit: dimensionless" << endl;
    file << "# xmin: 0" << endl;
    file << "# ymin: 0" << endl;
    file << "# zmin: 0" << endl;
    file << "# xmax: " << to_string(DELTA * N) << endl;
    file << "# ymax: " << to_string(DELTA * N) << endl;
    file << "# zmax: " << to_string(DELTA) << endl;
    file << "# valuedim: 3" << endl;
    file << "# valuelabels: m_x m_y m_z" << endl;
    file << "# Desc: Total simulation time: " << to_string(t) << endl;
    file << "# xbase: " << to_string(DELTA/2) << endl;
    file << "# ybase: " << to_string(DELTA/2) << endl;
    file << "# zbase: " << to_string(DELTA/2) << endl;
    file << "# xnodes: " << to_string(N) << endl;
    file << "# ynodes: " << to_string(N) << endl;
    file << "# znodes: 1" << endl;
    file << "# xstepsize: " << to_string(DELTA) << endl;
    file << "# ystepsize: " << to_string(DELTA) << endl;
    file << "# zstepsize: " << to_string(DELTA) << endl;
    file << "# End: Header" << endl;
    file << "# Begin: Data Text" << endl;

    // OVF uses the Fortran convention, i.e. x changes fastest
    for (int j=0; j<N; j++)
    {
        for (int i=0; i<N; i++)
        {
            file << m->data[m->getIndex(i, j, 0)] << " " << m->data[m->getIndex(i, j, 1)] << " " << m->data[m->getIndex(i, j, 2)] << endl;
        }
    }

    file << "# End: Data Text" << endl;
    file << "# End: Segment" << endl;

    file.close();
}
