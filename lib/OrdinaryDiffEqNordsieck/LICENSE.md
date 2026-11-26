The OrdinaryDiffEq.jl package is licensed under the MIT "Expat" License:

> Copyright (c) 2016-2020: ChrisRackauckas, Yingbo Ma, Julia Computing Inc, and
> other contributors:
> 
> https://github.com/SciML/OrdinaryDiffEq.jl/graphs/contributors
> 
> Permission is hereby granted, free of charge, to any person obtaining a copy
> of this software and associated documentation files (the "Software"), to deal
> in the Software without restriction, including without limitation the rights
> to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
> copies of the Software, and to permit persons to whom the Software is
> furnished to do so, subject to the following conditions:
> 
> The above copyright notice and this permission notice shall be included in all
> copies or substantial portions of the Software.
> 
> THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
> IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
> FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
> AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
> LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
> OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
> SOFTWARE.

## SUNDIALS

The functions `calc_coeff_bdf!` and `nordsieck_adjust_order_bdf!` in
`src/nordsieck_utils.jl` are translated from the functions `cvSetBDF`, `cvSetTqBDF`,
`cvIncreaseBDF` and `cvDecreaseBDF` in `src/cvode/cvode.c` of SUNDIALS v7.9.0
(https://github.com/LLNL/sundials/blob/v7.9.0/src/cvode/cvode.c). Those portions are
covered by the SUNDIALS license and notice below.

> BSD 3-Clause License
>
> Copyright (c) 2025-2026, Lawrence Livermore National Security, University of Maryland Baltimore County, and the SUNDIALS contributors.
> Copyright (c) 2013-2025, Lawrence Livermore National Security and Southern Methodist University.
> Copyright (c) 2002-2013, Lawrence Livermore National Security.
> All rights reserved.
>
> Redistribution and use in source and binary forms, with or without
> modification, are permitted provided that the following conditions are met:
>
> * Redistributions of source code must retain the above copyright notice, this
>   list of conditions and the following disclaimer.
>
> * Redistributions in binary form must reproduce the above copyright notice,
>   this list of conditions and the following disclaimer in the documentation
>   and/or other materials provided with the distribution.
>
> * Neither the name of the copyright holder nor the names of its
>   contributors may be used to endorse or promote products derived from
>   this software without specific prior written permission.
>
> THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
> AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
> IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
> DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
> FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
> DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
> SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
> CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
> OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
> OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

> This work was produced under the auspices of the U.S. Department of
> Energy by Lawrence Livermore National Laboratory under Contract
> DE-AC52-07NA27344.
>
> This work was prepared as an account of work sponsored by an agency of
> the United States Government. Neither the United States Government nor
> Lawrence Livermore National Security, LLC, nor any of their employees
> makes any warranty, expressed or implied, or assumes any legal liability
> or responsibility for the accuracy, completeness, or usefulness of any
> information, apparatus, product, or process disclosed, or represents that
> its use would not infringe privately owned rights.
>
> Reference herein to any specific commercial product, process, or service
> by trade name, trademark, manufacturer, or otherwise does not necessarily
> constitute or imply its endorsement, recommendation, or favoring by the
> United States Government or Lawrence Livermore National Security, LLC.
>
> The views and opinions of authors expressed herein do not necessarily
> state or reflect those of the United States Government or Lawrence
> Livermore National Security, LLC, and shall not be used for advertising
> or product endorsement purposes.
