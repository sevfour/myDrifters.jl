### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ f1498618-07d5-4892-9028-dabf977bff9b
begin
	using Pkg; Pkg.status()
		Pkg.develop(path="/Users/severinf/Scripts/git/Drifters.jl/")#to make sure it reloads when changed
	using Drifters, CairoMakie
	using Proj, MeshArrays, GeoJSON, PlutoUI, CSV, DataFrames
	lon0=-160.0; proj=Proj.Transformation(MA_preset=2,lon0=lon0)
	Pkg.status()
end

# ╔═╡ 39055364-6654-42dd-afad-9e67f286f054
md"""# Drifters.jl + Oscar Data Assimilative Model

- Oscar provides daily averaged surface currents are provided on a global 0.25 x 0.25 degree grid.
- Drifters.jl is used to compute trajectories of virtual floating parcels following Oscar flow fields.

!!! note
    See Appendix for more information on software and data.
"""

# ╔═╡ 1aed4830-43dc-4d98-bcc0-775738a477c1
md"""## Precompute Trajectories into a CSV File

!!! warning
    Define input path (raw oscar files) and precomputed output file (csv file)
"""

# ╔═╡ a8ddb075-73cf-4496-a8a5-4f824a5f80d6
begin
	input_path="oscar/data"
	input_files=Drifters.Oscar.list_files(input_path,2021)
	file_precomputed=joinpath(Drifters.datadeps.getdata("Oscar_2021_small"),"Drifters_Oscar_small.csv")
	println(input_files[1])
	println(input_files[end])
end

# ╔═╡ 3bb3a27c-a184-4257-8678-0eb33b980622
if !isempty(input_files)
	n_part=10000
	reset_rate=0.05
	nt=30
	I=Drifters.Oscar.main_loop(input_files=input_files,output_file=file_precomputed, n_part=n_part, reset_rate=reset_rate, nt=nt, do_save=false, verbose=true)
	I
end
# n_part = nb initial particles 
# reset_rate = reset per day
# nt= nb of days from the beginning we do calculation
# details is size (nt+1)*n_part

# ╔═╡ 2d690809-f527-4d9a-8a70-313100aa2927
begin
	df=I.🔴
	(t0,t1)=extrema(unique(df.t)[1:end]);
	CSV.write(file_precomputed,df) #store these data
	println(t0/86400)
	println(t1/86400)
end

# ╔═╡ 5e4f26a0-8954-44b7-a3cb-224197a2e0cb
md"""## Visualize Precomputed Data

- precomputed trajectory data are retrieved from a `csv` file.
- the file contains daily time series of positions and normalized velocities.
- virtual parcels were initially released the time varying Oscar flow fields.
- the parcel trajectories were then computed using Drifters.jl (code provided below).
"""

# ╔═╡ 1dc43b12-c49f-4ab5-a239-9188d8452659
# df=CSV.read(file_precomputed,DataFrame)

# ╔═╡ 08ba0bfd-b6b1-436a-8a7c-5e5f394fc465
begin
	times=sort(unique(df.t))
	tt=Observable(1) # particle position at that time seeded at the beginning of period (see csv file)
	nt0=10 # number of following positions to show for each particle (filaments)
	println(times[1:end]/86400) # in days
end

# ╔═╡ 7d345828-4078-4273-a62e-9e5f5354591d
begin
	🔴=@lift(filter(:t => x -> (x >= times[$tt])&&(x <=times[$tt+nt0-1] ), df))

	lon=@lift($🔴.lon)
	lat=@lift($🔴.lat)
	vel=@lift(86400*sqrt.($🔴.dxdt.^2+$🔴.dydt.^2))
	
	options=(plot_type=:Oscar_plot,proj=proj,lon0=-160,add_background=true,add_polygons=true,
			lon=lon,lat=lat,color=vel,colorrange=(0,2),colormap=:thermal,markersize=2)
	J=DriftersDataset( data=(df=df,), options=options)
	fig=plot(J)
end

# ╔═╡ 0845b7bf-2e77-4dde-a9a1-6d1ded585a01
begin
	file_output_mp4=tempname()*".mp4"
	record(fig, file_output_mp4, 1:nt-10, framerate = 25) do t
    	tt[]=t
	end
	print(file_output_mp4)
end

# ╔═╡ 24ecbc7b-b0f2-4ebe-9ac4-1a827254f225
md"""## Appendix

### Julia Packages
"""

# ╔═╡ 1f479c36-f021-464c-a279-0d62a1f33359
md"""### Software : Drifters.jl (v0.6.4)

Forget, G., (2021). IndividualDisplacements.jl: a Julia package to simulate and study particle displacements within the climate system. Journal of Open Source Software, 6(60), 2813, <https://doi.org/10.21105/joss.02813>

<https://github.com/JuliaClimate/Drifters.jl>

### Data : Oscar (v2.0)

ESR; Dohan, Kathleen. 2022. Ocean Surface Current Analyses Real-time (OSCAR) Surface Currents - Final 0.25 Degree (Version 2.0). Ver. 2.0. PO.DAAC, CA, USA. Dataset accessed [YYYY-MM-DD] at https://doi.org/10.5067/OSCAR-25F20
 
<https://podaac.jpl.nasa.gov/dataset/OSCAR_L4_OC_FINAL_V2.0>
  
#### Sample NetCDF Granule

<https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-protected/OSCAR_L4_OC_FINAL_V2.0/oscar_currents_final_20220504.nc>

#### Download Files

In shell : 

```
podaac-data-downloader -c OSCAR_L4_OC_FINAL_V2.0 -d ./data --start-date 2021-01-01T00:00:00Z --end-date 2021-02-01T00:00:00Z -e ""
```

#### Additional Information

- <https://github.com/podaac/data-subscriber>
- <https://podaac.jpl.nasa.gov/dataset/OSCAR_L4_OC_FINAL_V2.0#capability-modal-download>
"""

# ╔═╡ Cell order:
# ╟─39055364-6654-42dd-afad-9e67f286f054
# ╟─1aed4830-43dc-4d98-bcc0-775738a477c1
# ╠═a8ddb075-73cf-4496-a8a5-4f824a5f80d6
# ╠═3bb3a27c-a184-4257-8678-0eb33b980622
# ╠═2d690809-f527-4d9a-8a70-313100aa2927
# ╟─5e4f26a0-8954-44b7-a3cb-224197a2e0cb
# ╠═1dc43b12-c49f-4ab5-a239-9188d8452659
# ╠═08ba0bfd-b6b1-436a-8a7c-5e5f394fc465
# ╠═7d345828-4078-4273-a62e-9e5f5354591d
# ╠═0845b7bf-2e77-4dde-a9a1-6d1ded585a01
# ╟─24ecbc7b-b0f2-4ebe-9ac4-1a827254f225
# ╠═f1498618-07d5-4892-9028-dabf977bff9b
# ╟─1f479c36-f021-464c-a279-0d62a1f33359
