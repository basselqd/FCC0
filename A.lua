-- Simulation of 8x8x8 FCC array of nanoparticles
-- Maglua 307
--This file creat a macro lattice and calculate dipole field
size = mpi.get_size()
rank = mpi.get_rank()

n = 8 -- The lattice is nxnxn nanoparticle
-- the macro system (the point dipole interactions)
ssMacro = SpinSystem.new(n,n,n)
lr = LongRange3D.new(ssMacro)

lr:loadTensors("8x8x8_fcc.lua")

gx,gy,gz=n,n,n --dimentions of the lattice of nanosphere
dnumber=100  -- do dipole field update each dnumber steps


dfac= 3.69e-5  -- factor of multiplication of the dipole tensors (25 * Ub^2 * U0)/(4 kB Pi D^3) 
lr:setStrength(dfac)


if size ~= gx*gy*gz then
  error("This simulation can only be run with " .. gx*gy*gz .. " processors")
end


-- create mapping between rank (r) and macro lattice (nanoparticle position in the super lattice {i,j,k} ) 
mappingd = {}  
for k=0,gz-1 do
  for j=0,gy-1 do
    for i=0,gx-1 do
      r= 1+ (i)*n^0 +(j)*n^1 + (k)*n^2
      mappingd[r] = {i+1, j+1,k+1}
    end
  end
end



function fieldMacro() --calculate dipolar field at the macro level
  ssMacro:resetFields()
  lr:apply(ssMacro)
  ssMacro:sumFields()
end

ffm = io.open("dd.dat","a")

function writt() -- write the macro lattice info
  local en=0

  for rr=1,size do --calculate the dipole energy
    local coords = mappingd[rr]
    local mx,my,mz=ssMacro:spin(coords)
    local dx, dy, dz = ssMacro:field("Dipole", coord)
    en=en + mx*dx + my*dy +mz*dz
  end
  en = -en/(2*size)
  local tt=ss:time()
  local mx, my, mz,mt = ssMacro:netMoment(1/size)
  local lin=T.."\t"..tt.."\t"..mx.."\t"..my.."\t"..mz.."\t"..mt.."\t"..en.."\n"

  ffm:write(lin)
  ffm:flush()
end


allMoments = {} 

function applyFld() -- collect the net moment of each nanoparticle, calculate dipole field, apply it to each nanoparticle
  local mx, my, mz, mt = ss:netMoment() -- the moment of the current nanoparticle
  -- gather and send the moment of each nanoparticle to every processor
  allMoments = mpi.gather(1, {mx,my,mz})
  allMoments = mpi.bcast(1, allMoments)

  for ir=1,size do	--for all spheres
    local coords = mappingd[ir]
    ssMacro:setSpin(coords, allMoments[ir]) -- set the moment of the nanoparticle (ir) at the curresponding site in the macro array
  end
  fieldMacro() --calculate global  dipolar field
  local myCoord = mappingd[rank]
  local x1,y1,z1=mappingd[rank][1],mappingd[rank][2],mappingd[rank][3]
  local dx, dy, dz = ssMacro:field("Dipole",{x1,y1,z1})-- get the dipolar field on the sphere's site
  af:set({dx, dy, dz})-- apply/set dipole field on the current nanoparticle
end

-- set file for every nanoparticle to record info
cfn = string.format("dip."..mappingd[rank][1].."."..mappingd[rank][2].."."..mappingd[rank][3]..".dat", L)  
fd = io.open(cfn, "a")

function rayt() --write the individual nanoparticle info
  if rank==1 then -- to avoid all processors writing the same info
    writt()-- write the lattice info
  end
  local mx, my, mz, mt = ss:netMoment()
  local x1, y1, z1 = mappingd[rank][1],mappingd[rank][2],mappingd[rank][3] -- the position of the current nanoparicle
  local dx, dy, dz = ssMacro:field("Dipole",{x1,y1,z1}) --get the field at the position of the current nanoparticle
  local tt = ss:time()
	lin = T.."\t"..tt.."\t"..mx.."\t"..my.."\t"..mz.."\t"..mt.."\t"..dx.."\t"..dy.."\t"..dz.."\n"
	fd:write(lin)
	fd:flush()
end

--initial temperature
runTemp=70

if arg[1]=="load" then 
  T = arg[2] --set Temperature to argument [2]
  runTemp=T
  dol_t = 1 -- do load time for things that will be done only once at loading
  invmap, ss, ex, ani, af, mappingd, tTill2, next_func_call2, sstime, drun2, runTemp, rngsd=checkpointLoad("SSS"..T.."."..rank..".dat")
  dofile("DenseSphereL.lua")
else
  rngsd = rank
  dofile("DenseSphere.lua")
  dol_t = 0 -- do not load time for things that will be done only once at loading
end




function ryt_dtl(T,nm)  -- write configuration of current nanoparticle spins
  fc = io.open("core"..nm..T.."."..rank..".dat","a") -- for core spins
  fs = io.open("surf"..nm..T.."."..rank..".dat","a") -- for surface spins
  local ffx,ffy,ffz = 0,0,0
  local cct = 0
  for pos,ed in ss:extraDataIterator(false) do
    local cs = 0
    if ed.core then
      cs = "core"
    else
      cs = "surface"
    end
    local xi=pos[1]
    local sx,sy,sz,m=ss:spin({xi,1,1})
    local x,y,z=mapping[xi][1][1],mapping[xi][1][2],mapping[xi][1][3]
    local typ=ed.lbl
    local fx,fy,fz=ss:field("Exchange",{xi,1,1})
    local fix,fiy,fiz=ss:field("Anisotropy",{xi,1,1})
    ffx,ffy,ffz=fix+fx,fiy+fy,fiz+fz
    if  sx ~= 0 then
      ct = (fx*sx+fy*sy+fz*sz)/(fx^2+fy^2+fz^2)^(1/2)
    else
      ct = 0
    end

    cct = cct + ct
    local lin=x.."\t"..y.."\t"..z.."\t"..typ.."\t"..cs.."\t"..sx.."\t"..sy.."\t"..sz.."\t"..ffx.."\t"..ffy.."\t"..ffz.."\t"..cct.."\n"
    if ed.core then
      fc:write(lin)
      fc:flush()
    else
      fs:write(lin)
      fs:flush()
    end
  end
end



if rank == 1 then
  fdate=io.open("date.dat","a")
end

runt = 1000 -- default equilibration time


for di=runTemp,0,-5 do
  T = di
  if di<1 then
    T=0
  end

  if rank == 1 then
    fdate:write(T,"     ",ss:time(),"        ",os.date())
  end
  th:set(T)

  run(runt, 0.5)
  ryt_dtl(T,"c")

  if rank == 1 then
    fdate:write(T,"     ",ss:time(),"        ",os.date())
  end
end
