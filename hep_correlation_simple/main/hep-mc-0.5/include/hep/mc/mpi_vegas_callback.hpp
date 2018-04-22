#ifndef HEP_MC_MPI_VEGAS_CALLBACK_HPP
#define HEP_MC_MPI_VEGAS_CALLBACK_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2014-2015  Christopher Schwan
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "hep/mc/vegas_callback.hpp"
#include "hep/mc/vegas_result.hpp"

#include <functional>
#include <vector>

#include <mpi.h>

namespace hep
{

/// \addtogroup callbacks
/// @{

/// The default callback function. This function does nothing and always returns
/// `true`. It is the MPI equivalent of \ref vegas_default_callback.
///
/// \see mpi_vegas_callback
template <typename T>
inline bool mpi_vegas_default_callback(
	MPI_Comm,
	std::vector<vegas_result<T>> const&
) {
	return true;
}

/// Callback function that prints a detailed summary about every iteration
/// performed so far. This function always returns `true`. It is the equivalent
/// of \ref vegas_verbose_callback and only writes an output if it was called
/// from rank zero to avoid duplicated output.
///
/// \see vegas_callback
template <typename T>
inline bool mpi_vegas_verbose_callback(
	MPI_Comm communicator,
	std::vector<vegas_result<T>> const& results
) {
	int rank = -1;
	MPI_Comm_rank(communicator, &rank);

	if (rank == 0)
	{
		vegas_verbose_callback<T>(results);
	}

	return true;
}

/// The type of callback function that can be set by the user with
/// \ref mpi_vegas_callback.
template <typename T>
using mpi_vegas_callback_type =
	std::function<bool(MPI_Comm, std::vector<vegas_result<T>>)>;

/// Sets the vegas `callback` function and returns it. This function is called
/// after each iteration performed by \ref mpi_vegas. The default callback is
/// \ref mpi_vegas_default_callback which does nothing. The callback function
/// can e.g. be set to \ref mpi_vegas_verbose_callback which prints detailed
/// results after each iteration.
///
/// If this function is called without any argument, the current callback
/// function is returned.
template <typename T>
inline mpi_vegas_callback_type<T> mpi_vegas_callback(
	mpi_vegas_callback_type<T> callback = nullptr
) {
	static mpi_vegas_callback_type<T> object = mpi_vegas_default_callback<T>;

	if (callback != nullptr)
	{
		object = callback;
	}

	return object;
}

/// @}

}

#endif
