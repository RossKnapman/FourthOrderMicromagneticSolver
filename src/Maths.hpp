#pragma once

void crossProduct(float* (*a0), float* (*a1), float* (*a2), float* (*b0), float* (*b1), float* (*b2), float* (*out0), float* (*out1), float* (*out2))
{
    (*out0) = (*a1)*(*b2) - (*b1)*(*a2);
    (*out1) = (*a2)*(*b0) - (*a0)*(*b2);
    (*out2) = (*a0)*(*b1) - (*a1)*(*b0);
}
