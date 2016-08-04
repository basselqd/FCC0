-- This file set the micro system (nanoparticle) with stochastic
-- It is similar to DensSphere.lua with additions to load the system from a saved point
rank=mpi.get_rank()
rngsd=rank

dofile("maglua://RungeKutta.lua")
dofile("maglua://MakeMacro.lua")
dofile("SetupSphericalParticle.lua") -- load setup function
--dofile("DataCollection.lua") -- get magnetization functions

rng = Random.Isaac.new()
rng:setSeed(rngsd) 

-- setup system, objects, etc
nos,ss, ex, ani, regions, sinfo = makeSphericalParticle(L, ST, rng)

th  = Thermal.new(ss, rng)
af  = AppliedField.new(ss)
llg = LLG.Quaternion.new(ss)
llg:setThermalOnlyFirstTerm(false)
odt = ss:timeStep() -- get timestep

-- mapping here is to increase the efficiency of sLLG by mapping existing exchange only

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


if arg[1]=="load" then
   dol_t=1 -- do load time for things that will be done only once at loading
	invmap, ss, ex, ani, af, mappingd, tTill2, next_func_call2, sstime, drun2, runTemp, rngsd=checkpointLoad("SSS"..T.."."..rank..".dat")
   -- set new seed for rng and th
   rng=nil
   th=nil
   rng = Random.Isaac.new()
   rngsd=rngsd+1000
   rng:setSeed(rngsd) 
   T=runTemp
   ss:setTime(sstime)
   th  = Thermal.new(ss, rng)
end



ftime=io.open("time.dat","a")
lin=T.."\t"..rngsd.."\t"..ss:time().."\t"..os.date().."\t"..rank.."\n"
ftime:write(lin)



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


function run(dt, r)  --equilibration till tTill, and save info each r tu

   if dol_t==1 then
	ss:setTime(sstime)
	if arg[3]=="again" then
	   tTill=tonumber(arg[4])
	else
	   tTill=tTill2
	end
	next_func_call=next_func_call2
	drun=drun2
	dol_t=0
   else
	tTill=ss:time()+dt
	next_func_call=ss:time()+r
	drun=0
   end
   drun=0

   next_save_tim=ss:time()+40*r -- to make a save point every 40r
   sst = ss:time()
   local t=sst
   while ss:time() <  tTill do --for dt do
	DplFld() 
	step(ss)
	if ss:time() > next_func_call then -- each r time unites do ...
	   rayt()	
	   next_func_call = ss:time() + r
	   drun=0
	end

	if ss:time() > next_save_tim then -- each 40r time unites do ...
	   drun=0
	   checkpointSave("SSS"..T.."."..rank..".dat",invmap, ss, ex, ani, af, mappingd, tTill, next_func_call, ss:time(), drun, T,rngsd)
	   next_save_tim=ss:time()+40*r
	end
   end


   drun=0
   checkpointSave("SSS"..T..".f."..rank..".dat",invmap, ss, ex, ani, af, mappingd, tTill, next_func_call, ss:time(), drun, T,rngsd)

end




