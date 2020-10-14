!***  DESCRIPTION: this module define a number of interfaces to function 
!                  of cuRan. These routines can be used to generate random 
!                  number on GPU
!                  
!                 refer to the book "CUDA Fortran for Scientists and Engineers"
!                 __________________________________________________________
!                 ZHAI Lei, Der , 2017
!

module CudaRandomC2F_M
  use cudafor

  implicit none
  integer, public :: CURAND_RNG_PSEUDO_DEFAULT = 100 
  integer, public :: CURAND_RNG_PSEUDO_XORWOW  = 101 

  integer, public :: CURAND_RNG_QUASI_DEFAULT  = 200 
  integer, public :: CURAND_RNG_QUASI_SOBOL32  = 201
  integer, public :: CURAND_RNG_QUASI_SCRAMBLED_SOBOL32 = 202
  integer, public :: CURAND_RNG_QUASI_SOBOL64  = 203
  integer, public :: CURAND_RNG_QUASI_SCRAMBLED_SOBOL64 = 204

  interface curandCreateGenerator
     integer function curandCreateGenerator( &
          generator,rng_type) &
          bind(C,name='curandCreateGenerator') 
       use iso_c_binding
       integer(c_size_t):: generator
       integer(c_int),value:: rng_type
     end function curandCreateGenerator
  end interface curandCreateGenerator
  
  interface curandSetPseudoRandomGeneratorSeed
     integer function curandSetPseudoRandomGeneratorSeed( &
          generator,seed) &
          bind(C,name='curandSetPseudoRandomGeneratorSeed')
       use iso_c_binding
       integer(c_size_t), value:: generator
       integer(c_long_long),value:: seed
     end function curandSetPseudoRandomGeneratorSeed
  end interface curandSetPseudoRandomGeneratorSeed
  
  interface curandGenerateUniform
     integer function curandGenerateUniform( &
          generator, odata, numele) &
          bind(C,name='curandGenerateUniform')
       use iso_c_binding
       integer(c_size_t),value:: generator
       !pgi$ ignore_tr odata
       real(kind=4), device:: odata(*)
       integer,value:: numele
     end function curandGenerateUniform
     
     integer function curandGenerateUniformDouble(&
          generator, odata, numele) &
          bind(C,name='curandGenerateUniformDouble')
       use iso_c_binding
       integer(c_size_t),value:: generator
       !pgi$ ignore_tr odata
       real(kind=8), device:: odata(*)
       integer,value:: numele
     end function curandGenerateUniformDouble
  end interface curandGenerateUniform
  
  interface curandGenerateNormal
     integer function curandGenerateNormal( &
          generator, odata, numele, mean,stddev) &
          bind(C,name='curandGenerateNormal')
       use iso_c_binding
       integer(c_size_t),value:: generator
       !pgi$ ignore_tr odata
       real(kind=4), device:: odata(*)
       integer,value:: numele
       real(c_float), value:: mean,stddev
     end function curandGenerateNormal
     
     integer function curandGenerateNormalDouble( &
          generator, odata, numele,mean, stddev) &
          bind(C,name='curandGenerateNormalDouble')
       use iso_c_binding
       integer(c_size_t),value:: generator
       !pgi$ ignore_tr odata
       real(kind=8), device:: odata(*)
       integer,value:: numele
       real(c_double), value:: mean,stddev
     end function curandGenerateNormalDouble
  end interface curandGenerateNormal
  
  interface curandDestroyGenerator
     integer function curandDestroyGenerator(generator) &
          bind(C,name='curandDestroyGenerator')
       use iso_c_binding
       integer(c_size_t),value:: generator
     end function curandDestroyGenerator
  end interface curandDestroyGenerator
  
end module CudaRandomC2F_M