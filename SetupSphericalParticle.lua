-- This script/function sets up a spherical particle and its exchange and anisotropy
--         requires: maglua-r307 or newer
-- 

-- Input:
--  L = radius (unit cells)
-- 
--rng = random number generator
L   = 9

function makeSphericalParticle(L, ST, rng)
	local RS = 4*L
	local RC = RS*(6.75/7.5)


	local kcore = 0.0 * 2 -- core anisotropy
	local ksurf = 10      -- surface anisotropy
-- core-core exchange   25 comes from normalizing the spin 5uB*5uB
	local rjaacc=-21*25
	local rjabcc=-28.1*25
	local rjbbcc=-8.6*25

	local  acs=0.05/2  -- exchange ratio of  (core-surface)/(core-core)
	local  ass=0.05/2 -- exchange ratio of  (surface-surface)/(core-core)
-- core-surface exchange
	local rjaacs=rjaacc*acs
	local rjabcs=rjabcc*acs
	local rjbbcs=rjbbcc*acs
-- surface-surface exchange
	local rjaass=rjaacc*ass
	local rjabss=rjabcc*ass
	local rjbbss=rjbbcc*ass

	-- used for initial state stats. matching MC code
	local nmax, nc, ns, nv = 0, 0, 0, 0
	local types = {"A1", "A2", "B1", "B2", "B3", "B4"}
	-- this will hold all site information during initialization
	local sites = {} 

	-- given a coordinate, classify site and add to temporary sites table
	local function addSite(x, y, z, lbl, occ)
		if x^2 + y^2 + z^2 <= RS^2 then
			local s={
					x=x, y=y, z=z,
					lbl =lbl,
					core=(x^2 + y^2 + z^2) <= RC^2,
					vac =rng:rand() > occ,
					radial = {x,y,z}
				}
			table.insert(sites, s)
		
			nmax = nmax + 1
			if s.core then
				nc = nc + 1
			else
				ns = ns + 1
			end
			if s.vac then
				nv = nv + 1
			end
		end
	end

	-- populating the sparse cubic system with sites as per
	-- the pattern in Byron's MC code
	for i=-L,L do
		for j=-L,L do
			for k=-L,L do
				-- A1 sites (0,0,0)+fcc
				local ax1 = 4*(i+k)
				local ay1 = 4*(i+j)
				local az1 = 4*(j+k)
				addSite(ax1, ay1, az1, "A1", 1.0)
				
				-- A1 sites (1/4,1/4,1/4)+fcc
				local ax2 = ax1+2
				local ay2 = ay1+2
				local az2 = az1+2
				addSite(ax2, ay2, az2, "A2", 1.0)
				
				-- B1 sites (1/8,5/8,1/8)+fcc
				local bx1 = ax1 + 1
				local by1 = ay1 + 1
				local bz1 = az1 - 3
				addSite(bx1, by1, bz1, "B1", 5/6)
				
				-- B2 sites (3/8,5/8,3/8)+fcc
				local bx2 = ax1 + 3
				local by2 = ay1 + 1
				local bz2 = az1 - 1
				addSite(bx2, by2, bz2, "B2", 5/6)
				
				-- B3 sites (3/8,7/8,1/8)+fcc
				local bx3 = ax1 - 1
				local by3 = ay1 + 3
				local bz3 = az1 + 1
				addSite(bx3, by3, bz3, "B3", 5/6)
				
				-- B4 sites (1/8,7/8,3/8)+fcc
				local bx4 = ax1 + 1
				local by4 = ay1 + 3
				local bz4 = az1 - 1
				addSite(bx4, by4, bz4, "B4", 5/6)
			end
		end
	end

	local info = "#sites="..nmax.." #core="..nc.." #surface="..ns.." #vacancies="..nv

	-- These are the offsets to the neighbours for each site type
	local nn = {}
	nn["A1"] = {{ 2, 2, 2}, { 2,-2,-2}, {-2, 2,-2}, {-2,-2, 2}, -- a2
				{ 1, 1,-3}, { 1,-3, 1}, {-3, 1, 1}, -- b1
				{-1, 1, 3}, {-1,-3,-1}, { 3, 1,-1}, -- b2
				{-1,-1,-3}, {-1, 3, 1}, { 3,-1, 1}, -- b3
				{ 1,-1, 3}, { 1, 3,-1}, {-3,-1,-1}} -- b4

	nn["A2"] = {{-2,-2,-2}, {-2, 2, 2}, { 2,-2, 2}, { 2, 2,-2}, -- a1
				{-1,-1, 3}, {-1, 3,-1}, { 3,-1,-1}, -- b1
				{ 1,-1,-3}, { 1, 3, 1}, {-3,-1, 1}, -- b2
				{ 1, 1, 3}, { 1,-3,-1}, {-3, 1,-1}, -- b3
				{-1, 1,-3}, {-1,-3, 1}, { 3, 1, 1}} -- b4

	nn["B1"] = {{ 2, 0, 2}, {-2, 0,-2}, { 2, 2, 0}, -- b2,b3,b4
				{-2,-2, 0}, { 0, 2, 2}, { 0,-2,-2}, 
				{-1,-1, 3}, {-1, 3,-1}, { 3,-1,-1}, -- a1,a2
				{ 1, 1,-3}, { 1,-3, 1}, {-3, 1, 1}}

	nn["B2"] = {{ 2, 0, 2}, {-2, 0,-2}, {-2, 2, 0}, -- b1,b3,b4
				{ 2,-2, 0}, { 0, 2,-2}, { 0,-2, 2},
				{ 1,-1,-3}, { 1, 3, 1}, {-3,-1, 1}, -- a1,a2
				{-1, 1, 3}, {-1,-3,-1}, { 3, 1,-1}}

	nn["B3"] = {{ 2, 0,-2}, {-2, 0, 2}, { 2, 2, 0}, -- b1,b2,b4
				{-2,-2, 0}, { 0, 2,-2}, { 0,-2, 2},
				{ 1, 1, 3}, { 1,-3,-1}, {-3, 1,-1}, -- a1,a2
				{-1,-1,-3}, {-1, 3, 1}, { 3,-1, 1}}

	nn["B4"] = {{ 2, 0,-2}, {-2, 0, 2}, { 2,-2, 0}, -- b1,b2,b3
				{-2, 2, 0}, { 0, 2, 2}, { 0,-2,-2},
				{-1, 1,-3}, {-1,-3, 1}, { 3, 1, 1}, -- a1,a2
				{ 1,-1, 3}, { 1, 3,-1}, {-3,-1,-1}}

	-- shift all sites so they have positive coordinates
	local minx, miny, minz = 1, 1, 1
	local maxx, maxy, maxz = 1, 1, 1
	for k,v in pairs(sites) do
		minx = math.min(minx, v.x)                                                   ------------
		miny = math.min(miny, v.y) 
		minz = math.min(minz, v.z) 
		
		maxx = math.max(maxx, v.x) 
		maxy = math.max(maxy, v.y) 
		maxz = math.max(maxz, v.z) 
	end

	for k,v in pairs(sites) do
		sites[k].x = sites[k].x - (minx - 1)
		sites[k].y = sites[k].y - (miny - 1)
		sites[k].z = sites[k].z - (minz - 1)
	end

	maxx = maxx - minx + 1
	maxy = maxy - miny + 1
	maxz = maxz - minz + 1

	-- Create data and operator objects

	local ss  = SpinSystem.new(maxx, maxy, maxz)
	local ex  = Exchange.new(ss)
	local ani = Anisotropy.new(ss)

	-- now we will put sites in extraData and clear the table
	for k,s in pairs(sites) do
		local pos = {s.x, s.y, s.z}
		ss:setExtraData(pos, s)
	end
	sites = nil --the garbage collector will recover this big chunk of memory

	-- calculating neighbour sites for all sites (via extra data iterator)
	for pos,s in ss:extraDataIterator(false) do
		local neighbours = {}
		local x, y, z = s.x, s.y, s.z
		for k,v in pairs(nn[s.lbl]) do
			table.insert(neighbours, {x+v[1], y+v[2], z+v[3]}) 
		end
		s.neighbours = neighbours
		ss:setExtraData(pos, s)
	end

	-- lookup exchange strengths based on site types and core/surf
	function exStr(s1, s2)
		local t1 = string.sub(s1.lbl,1,1) -- first letter of site label
		local t2 = string.sub(s2.lbl,1,1)
		local c1 = s1.core
		local c2 = s2.core
		
		if t1 == "A" and t2 == "A" then
			if c1 ~= c2 then return rjaacs end --c/s
			if c1       then return rjaacc end --c/c
			return rjaass                      --s/s
		end

		if t1 == "B" and t2 == "B" then
			if c1 ~= c2 then return rjbbcs end
			if c1       then return rjbbcc end
			return rjbbss
		end

		-- A-B
		if c1 ~= c2 then return rjabcs end 
		if c1       then return rjabcc end
		return rjabss
	end

	-- does site exist and is not a vacancy?
	function real_site(x,y,z)
		if type(x) == "table" then -- this function can accept a table and cast as 3 numbers
			return real_site(x[1] or 1, x[2] or 1, x[3] or 1)
		end
		if x < 1 or y < 1 or z < 1 then
			return false
		end
		
		if x > maxx or y > maxy or z > maxz then
			return false
		end

		if not ss:extraData(x,y,z) then
			return false
		end
		
		return not ss:extraData(x,y,z).vac
	end

	-- setup system: initial orientation, exchage & anisotropy
	for z=1,maxz do
		for y=1,maxy do
			for x=1,maxx do
				if real_site(x,y,z) then
					local s1 = ss:extraData(x,y,z)

					-- initial orientation
					--  site {x,y,z} will point in a random 
					--  direction with unit magnetization
					ss:setSpin({x,y,z}, 
							{rng:normal(),rng:normal(),rng:normal()}, 1)
					
					-- setup exchange interaction
					for k,v in pairs(s1.neighbours) do
						local a, b, c = v[1], v[2], v[3]
						if real_site(a,b,c) then
							local s2 = ss:extraData(a,b,c)
							ex:addPath({x,y,z}, {a,b,c}, exStr(s1, s2))
						end
					end
			
					-- setup anisotropy
					if s1.core then
						ani:add({x,y,z}, {0,0,1}, kcore)
					else
						ani:add({x,y,z}, s1.radial, ksurf)
					end
				end
			end
		end
	end

	ss:setTimeStep(0.0002)
	ss:setAlpha(0.5)



	nos=nmax-nv
	return nos,ss, ex, ani, info
end
	

