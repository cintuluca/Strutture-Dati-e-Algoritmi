using ITensors
let
  # Create 100 spin-one indices
  N = 100
  sites = siteinds("S=1",N)

  # Input operator terms which define
  # a Hamiltonian matrix, and convert
  # these terms to an MPO tensor network
  # (here we make the 1D Heisenberg model)
  os = OpSum()
  for j=1:N-1
    os += "Sz",j,"Sz",j+1
    os += 0.5,"S+",j,"S-",j+1
    os += 0.5,"S-",j,"S+",j+1
  end
  H = MPO(os,sites)

  # Create an initial random matrix product state
  psi0 = randomMPS(sites)

  # Plan to do 5 passes or 'sweeps' of DMRG,
  # setting maximum MPS internal dimensions
  # for each sweep and maximum truncation cutoff
  # used when adapting internal dimensions:
  nsweeps = 5
  maxdim = [10,20,100,100,200]
  cutoff = 1E-10

  # Run the DMRG algorithm, returning energy
  # (dominant eigenvalue) and optimized MPS
  energy, psi = dmrg(H,psi0; nsweeps, maxdim, cutoff)
  println("Final energy = $energy")

  nothing
end
