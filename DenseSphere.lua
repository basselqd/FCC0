-- This script set the micro system (nanoparticle) with sLLG 
-- This script runs only once.  
rank=mpi.get_rank()
rngsd=rank+1000000 -- seed for random number generator

dofile("maglua://RungeKutta.lua")
dofile("maglua://MakeMacro.lua")
dofile("SetupSphericalParticle.lua") -- load setup function

rng = Random.Isaac.new()
rng:setSeed(rngsd)

-- setup system, objects, etc
nos,ss, ex, ani, regions, sinfo = makeSphericalParticle(L, ST, rng)

th  = Thermal.new(ss, rng)
af  = AppliedField.new(ss)
llg = LLG.Quaternion.new(ss)
llg:setThermalOnlyFirstTerm(false)
odt = ss:timeStep() --original timestep



-- mapping here is to increase the efficiency of sLLG 
mapping = {}

for k=1,ss:nz() do
	for j=1,ss:ny() do
		for i=1,ss:nx() do
			local ed = ss:extraData(i,j,k)
			if ed then --if the site exists
					table.insert(mapping, {{i,j,k}}) --  mapping for all
			end
		end
	end
end

-- combine, make new data/operators as per mapping
invmap, ss, ex, ani, af = MakeMacro(mapping, ss, ex, ani, af)


ftime=io.open("time.dat","a")

function calculateDeterministicField(spinsystem)
	spinsystem:resetFields()
	ani:apply(spinsystem)
	ex:apply(spinsystem)
	af:apply(spinsystem)
	spinsystem:sumFields()
end


function dynamics(spinsystem)
	local t = spinsystem:time()
	th:set(T)
end

-- sLLG step
step = make_rk_step_function(ss, "RK4", calculateDeterministicField, dynamics, llg, th)


function DplFld() --calculate and set the dipole field
	if (math.abs(math.fmod (drun,dnumber)) < 1 ) then --update the dipole field each dnumber (100) steps	
	applyFld()
	end
	drun=drun+1 -- counter for updating the dipole field
end

function run(dt, r) --equilibration for dt time units (tu) and save info each r tu

   tTill=ss:time()+dt
   next_func_call=ss:time()+r 
   drun=0 

	 sst = ss:time()
	local t=sst
	while ss:time() <  tTill do -- equilibrate till tTill
		DplFld() 
		step(ss)
		if ss:time() > next_func_call then -- each r time unites do ...
		   rayt()			   -- write individual nanoparticle info	

		   next_func_call = ss:time() + r
		   drun=0 -- update dipole field for the next sLLG step  (this is just for consistency in case of loading files later on)
		   -- save point
		   checkpointSave("SSS"..T.."."..rank..".dat",invmap, ss, ex, ani, af, mappingd, tTill, next_func_call, ss:time(), drun, T,rngsd)

		end
	end

   drun=0
   -- save the nanoparticle info and state to a file
   checkpointSave("SSS"..T..".f."..rank..".dat",invmap, ss, ex, ani, af, mappingd, tTill, next_func_call, ss:time(), drun, T,rngsd)
end


