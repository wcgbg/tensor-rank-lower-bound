#pragma once

// A 3-tensor with compile-time dimensions over 𝔽_q (q = P^M), with one
// GF<P, M> coefficient per cell.

#include <array>
#include <cstddef>

#include "core/gf.h"

template <int P, int M, std::size_t NA, std::size_t NB, std::size_t NC>
using Tensor = std::array<std::array<std::array<GF<P, M>, NC>, NB>, NA>;
