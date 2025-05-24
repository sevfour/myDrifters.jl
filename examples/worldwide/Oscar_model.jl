### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ f1498618-07d5-4892-9028-dabf977bff9b
begin
	using Pkg; Pkg.status()
	Pkg.develop(path="/Users/severinf/Scripts/git/Drifters.jl/")#to make sure it reloads when changed?
	using CairoMakie, Drifters
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

# ╔═╡ e278df98-196b-4f05-a4e1-010b287d221f
TableOfContents()

# ╔═╡ 5e4f26a0-8954-44b7-a3cb-224197a2e0cb
md"""## Visualize Precomputed Data

- 30 days of precomputed trajectory data are retrieved from a `csv` file.
- the file contains daily time series of positions and normalized velocities.
- 50000 virtual parcels were initially released the time varying Oscar flow fields.
- the parcel trajectories were then computed using Drifters.jl (code provided below).
"""

# ╔═╡ d6bd64ce-ebd3-11ef-3ea6-c77f6094d0a6
file_precomputed=joinpath(Drifters.datadeps.getdata("Oscar_2021_small"),"Drifters_Oscar_small.csv")

# ╔═╡ 1dc43b12-c49f-4ab5-a239-9188d8452659
df=CSV.read(file_precomputed,DataFrame)

# ╔═╡ 7d345828-4078-4273-a62e-9e5f5354591d
begin #let
	times=sort(unique(df.t))
	t0=Observable(1)
	nt0=2
	🔴=@lift(filter(:t => x -> (x >= times[$t0])&&(x <=times[$t0+nt0-1] ), df))

	lon=@lift($🔴.lon)
	lat=@lift($🔴.lat)
	vel=@lift(86400*sqrt.($🔴.dxdt.^2+$🔴.dydt.^2))
	
	options=(plot_type=:Oscar_plot,proj=proj,lon0=-160,add_background=true,add_polygons=true,
			lon=lon,lat=lat,color=vel,colorrange=(0,2),colormap=:thermal,markersize=2)
	J=DriftersDataset( data=(df=df,), options=options)
	fig=plot(J)
end

# ╔═╡ c187acd0-50a5-4360-ad13-73d658ebe87f
# Animation
begin
	nt=30
	file_output_mp4=tempname()*".mp4"
	record(fig, file_output_mp4, 1:nt-10, framerate = 25) do t
    	t0[]=t
	end
	print(file_output_mp4)
end

# ╔═╡ 1aed4830-43dc-4d98-bcc0-775738a477c1
md"""## Compute New Trajectories

!!! warning
    Recomputing trajectories requires Oscar data to have been downloaded to the `input_path` folder.
"""

# ╔═╡ a8ddb075-73cf-4496-a8a5-4f824a5f80d6
begin
	input_path="oscar/data"
	input_files=Drifters.Oscar.list_files(input_path,2021)
end

# ╔═╡ dd7f5d46-6479-4612-a6a4-2fedca50b1a8
if !isempty(input_files)
	n_part=10000
	reset_rate=0.05
	# nt=30
	I=Drifters.Oscar.main_loop(input_files=input_files,n_part=n_part, reset_rate=reset_rate, nt=nt, do_save=false, verbose=true)

	# options=(plot_type=:Oscar_plot,proj=proj,lon0=-160,add_background=true,add_polygons=true,markersize=2)
	Jnew=DriftersDataset( data=(df=I.🔴,), options=options);
	"all set"
end

# ╔═╡ 2a360d87-29c9-44f8-a467-f2648ce53948
isempty(input_files) ? nothing : plot(Jnew)

# ╔═╡ 24ecbc7b-b0f2-4ebe-9ac4-1a827254f225
md"""## Appendix

### Julia Packages
"""

# ╔═╡ 4a122475-11ab-4d1e-9ce6-c12148339430
isdir(MeshArrays.GRID_LLC90) ? "all set for plotting" : "missing background data for plotting"

# ╔═╡ 1f479c36-f021-464c-a279-0d62a1f33359
md"""### Software : Drifters.jl (v0.6.4)

Forget, G., (2021). IndividualDisplacements.jl: a Julia package to simulate and study particle displacements within the climate system. Journal of Open Source Software, 6(60), 2813, <https://doi.org/10.21105/joss.02813>

<https://github.com/JuliaClimate/Drifters.jl>

### Data : Oscar (v2.0)

To compute trajectories you'll want to : 

- download Oscar data from <https://podaac.jpl.nasa.gov/dataset/OSCAR_L4_OC_FINAL_V2.0>.
- put them in a subfolder called `oscar/data/` that should look as shown below.

```
OSCAR_L4_OC_FINAL_V2.0.citation.txt
oscar_currents_final_19930101.nc
oscar_currents_final_19930101.nc.md5
...
```

!!! note
    Data set citation : ESR; Dohan, Kathleen. 2022. Ocean Surface Current Analyses Real-time (OSCAR) Surface Currents - Final 0.25 Degree (Version 2.0). Ver. 2.0. PO.DAAC, CA, USA. Dataset accessed [YYYY-MM-DD] at https://doi.org/10.5067/OSCAR-25F20

#### More on Oscar

- Sample data granule : [this link](https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-protected/OSCAR_L4_OC_FINAL_V2.0/oscar_currents_final_20220504.nc)
- How-to [download](<https://podaac.jpl.nasa.gov/dataset/OSCAR_L4_OC_FINAL_V2.0#capability-modal-download) files : in a linux terminal, for example, use [podaac-data-downloader](https://github.com/podaac/data-subscriber)

```
podaac-data-downloader -c OSCAR_L4_OC_FINAL_V2.0 -d ./data --start-date 2021-01-01T00:00:00Z --end-date 2021-02-01T00:00:00Z -e ""
```
"""

# ╔═╡ Cell order:
# ╟─39055364-6654-42dd-afad-9e67f286f054
# ╟─e278df98-196b-4f05-a4e1-010b287d221f
# ╟─5e4f26a0-8954-44b7-a3cb-224197a2e0cb
# ╠═d6bd64ce-ebd3-11ef-3ea6-c77f6094d0a6
# ╠═1dc43b12-c49f-4ab5-a239-9188d8452659
# ╠═7d345828-4078-4273-a62e-9e5f5354591d
# ╠═c187acd0-50a5-4360-ad13-73d658ebe87f
# ╟─1aed4830-43dc-4d98-bcc0-775738a477c1
# ╠═a8ddb075-73cf-4496-a8a5-4f824a5f80d6
# ╠═dd7f5d46-6479-4612-a6a4-2fedca50b1a8
# ╠═2a360d87-29c9-44f8-a467-f2648ce53948
# ╟─24ecbc7b-b0f2-4ebe-9ac4-1a827254f225
# ╠═f1498618-07d5-4892-9028-dabf977bff9b
# ╟─4a122475-11ab-4d1e-9ce6-c12148339430
# ╟─1f479c36-f021-464c-a279-0d62a1f33359
