Overhead ╎ [+additional indent] Count File:Line; Function
=========================================================
     ╎50     @Base/array.jl:376; _copyto_impl!(dest::Matrix{Float64}, doffs::In…
     ╎ 50     @Base/array.jl:337; unsafe_copyto!
   50╎  50     @Base/cmem.jl:26; memmove
     ╎11     @Base/array.jl:2041; vcat(::Vector{Float64}, ::Vector{Float64})
     ╎ 11     @Base/array.jl:337; unsafe_copyto!
   11╎  11     @Base/cmem.jl:26; memmove
     ╎132306 @Base/client.jl:552; _start()
     ╎ 132306 @Base/client.jl:333; exec_options(opts::Base.JLOptions)
     ╎  132306 @Base/client.jl:416; run_main_repl(interactive::Bool, quiet::Boo…
     ╎   132306 @Base/essentials.jl:889; invokelatest
     ╎    132306 @Base/essentials.jl:892; #invokelatest#2
     ╎     132306 @Base/client.jl:432; (::Base.var"#1013#1015"{Bool, Bool, Bool…
     ╎    ╎ 132306 …0/REPL/src/REPL.jl:375; run_repl(repl::REPL.AbstractREPL, c…
     ╎    ╎  132306 …/REPL/src/REPL.jl:389; run_repl(repl::REPL.AbstractREPL, c…
     ╎    ╎   132306 …/REPL/src/REPL.jl:228; kwcall(::NamedTuple, ::typeof(REPL…
     ╎    ╎    132306 …/REPL/src/REPL.jl:231; start_repl_backend(backend::REPL.…
     ╎    ╎     132306 …REPL/src/REPL.jl:246; repl_backend_loop(backend::REPL.R…
     ╎    ╎    ╎ 132306 …REPL/src/REPL.jl:150; eval_user_input(ast::Any, backen…
     ╎    ╎    ╎  132306 @Base/boot.jl:385; eval
     ╎    ╎    ╎   132306 @Base/client.jl:489; include(fname::String)
     ╎    ╎    ╎    132306 …ase/loading.jl:2136; _include(mapexpr::Function, mo…
     ╎    ╎    ╎     132306 …ase/loading.jl:2076; include_string(mapexpr::typeo…
     ╎    ╎    ╎    ╎ 132306 @Base/boot.jl:385; eval
     ╎    ╎    ╎    ╎  132306 …nce_eclipse.jl:10; kwcall(::@NamedTuple{verbose:…
     ╎    ╎    ╎    ╎   132306 …ce_eclipse.jl:18; synthesize_spectra_eclipse(sp…
     ╎    ╎    ╎    ╎    3      …ce_eclipse.jl:35; synth_Eclipse_cpu(spec::Spec…
     ╎    ╎    ╎    ╎     3      …aceEclipse.jl:25; SynthWorkspaceEclipse
     ╎    ╎    ╎    ╎    ╎ 1      …ceEclipse.jl:40; GRASS.SynthWorkspaceEclipse…
     ╎    ╎    ╎    ╎    ╎  1      …se/array.jl:631; zeros(::Int64, ::Int64)
     ╎    ╎    ╎    ╎    ╎   1      …se/array.jl:633; zeros
     ╎    ╎    ╎    ╎    ╎    1      …se/array.jl:637; zeros
     ╎    ╎    ╎    ╎    ╎     1      …se/array.jl:395; fill!
     ╎    ╎    ╎    ╎    ╎    ╎ 1      …se/array.jl:1021; setindex!
     ╎    ╎    ╎    ╎    ╎ 1      …ceEclipse.jl:43; GRASS.SynthWorkspaceEclipse…
     ╎    ╎    ╎    ╎    ╎  1      …se/array.jl:631; zeros(::Int64, ::Int64)
     ╎    ╎    ╎    ╎    ╎   1      …se/array.jl:633; zeros
     ╎    ╎    ╎    ╎    ╎    1      …se/array.jl:637; zeros
     ╎    ╎    ╎    ╎    ╎     1      …se/array.jl:395; fill!
    1╎    ╎    ╎    ╎    ╎    ╎ 1      …se/array.jl:1021; setindex!
     ╎    ╎    ╎    ╎    ╎ 1      …ceEclipse.jl:46; GRASS.SynthWorkspaceEclipse…
     ╎    ╎    ╎    ╎    ╎  1      …rraymath.jl:356; repeat(::Vector{Tuple{Symb…
     ╎    ╎    ╎    ╎    ╎   1      …rraymath.jl:398; repeat
     ╎    ╎    ╎    ╎    ╎    1      …rraymath.jl:401; repeat(arr::Vector{Tuple…
     ╎    ╎    ╎    ╎    ╎     1      …rraymath.jl:459; repeat_inner_outer
     ╎    ╎    ╎    ╎    ╎    ╎ 1      …rraymath.jl:465; repeat_outer(a::Matrix…
     ╎    ╎    ╎    ╎    ╎    ╎  1      …actarray.jl:832; similar
     ╎    ╎    ╎    ╎    ╎    ╎   1      …se/array.jl:420; similar
     ╎    ╎    ╎    ╎    ╎    ╎    1      …ase/boot.jl:487; Array
    1╎    ╎    ╎    ╎    ╎    ╎     1      …ase/boot.jl:479; Array
     ╎    ╎    ╎    ╎    1      …ce_eclipse.jl:36; synth_Eclipse_cpu(spec::Spec…
     ╎    ╎    ╎    ╎     1      …aceEclipse.jl:24; GRASS.GeoWorkspaceEclipse(d…
     ╎    ╎    ╎    ╎    ╎ 1      …ase/array.jl:631; zeros(::Int64, ::Int64)
     ╎    ╎    ╎    ╎    ╎  1      …se/array.jl:633; zeros
     ╎    ╎    ╎    ╎    ╎   1      …se/array.jl:637; zeros
     ╎    ╎    ╎    ╎    ╎    1      …se/array.jl:395; fill!
    1╎    ╎    ╎    ╎    ╎     1      …se/array.jl:1021; setindex!
     ╎    ╎    ╎    ╎    1      …ce_eclipse.jl:39; synth_Eclipse_cpu(spec::Spec…
     ╎    ╎    ╎    ╎     1      …ase/array.jl:636; zeros(::Type{Int64}, dims::…
     ╎    ╎    ╎    ╎    ╎ 1      @Base/boot.jl:487; Array
    1╎    ╎    ╎    ╎    ╎  1      …ase/boot.jl:479; Array
     ╎    ╎    ╎    ╎    217    …ce_eclipse.jl:54; synth_Eclipse_cpu(spec::Spec…
     ╎    ╎    ╎    ╎     217    …/SolarData.jl:22; kwcall(::@NamedTuple{fname:…
     ╎    ╎    ╎    ╎    ╎ 217    …SolarData.jl:40; SolarData(; fname::String, …
     ╎    ╎    ╎    ╎    ╎  217    …SolarData.jl:46; SolarData
     ╎    ╎    ╎    ╎    ╎   216    …olarData.jl:62; SolarData(fname::String; r…
     ╎    ╎    ╎    ╎    ╎    216    …src/file.jl:94; h5open
     ╎    ╎    ╎    ╎    ╎     1      …src/file.jl:95; #h5open#16
     ╎    ╎    ╎    ╎    ╎    ╎ 1      …src/file.jl:63; h5open(filename::String…
     ╎    ╎    ╎    ╎    ╎    ╎  1      …src/file.jl:75; h5open(filename::Strin…
     ╎    ╎    ╎    ╎    ╎    ╎   1      …src/file.jl:20; h5open
     ╎    ╎    ╎    ╎    ╎    ╎    1      …src/file.jl:48; h5open(filename::Str…
     ╎    ╎    ╎    ╎    ╎    ╎     1      …src/file.jl:188; ishdf5(name::Strin…
    1╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …unctions.jl:1548; h5f_is_hdf5(path…
     ╎    ╎    ╎    ╎    ╎     215    …src/file.jl:96; #h5open#16
     ╎    ╎    ╎    ╎    ╎    ╎ 215    …ase/task.jl:297; task_local_storage(bod…
     ╎    ╎    ╎    ╎    ╎    ╎  213    …src/file.jl:101; (::HDF5.var"#17#18"{H…
     ╎    ╎    ╎    ╎    ╎    ╎   4      …olarData.jl:70; (::GRASS.var"#16#25"{…
     ╎    ╎    ╎    ╎    ╎    ╎    4      …ighlevel.jl:74; getindex
     ╎    ╎    ╎    ╎    ╎    ╎     1      …ighlevel.jl:75; getindex(parent::HD…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …c/groups.jl:128; haskey(parent::HD…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …c/groups.jl:132; haskey(parent::H…
    1╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …unctions.jl:2076; h5l_exists(loc…
     ╎    ╎    ╎    ╎    ╎    ╎     3      …ighlevel.jl:77; getindex(parent::HD…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 3      …/objects.jl:37; open_object
    3╎    ╎    ╎    ╎    ╎    ╎    ╎  3      …unctions.jl:2562; h5o_open(loc_id…
     ╎    ╎    ╎    ╎    ╎    ╎   1      …olarData.jl:91; (::GRASS.var"#16#25"{…
     ╎    ╎    ╎    ╎    ╎    ╎    1      …c/groups.jl:146; keys(x::HDF5.Group)
     ╎    ╎    ╎    ╎    ╎    ╎     1      …/helpers.jl:487; h5l_iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …/helpers.jl:491; h5l_iterate(f::An…
    1╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …unctions.jl:2143; h5l_iterate(gro…
     ╎    ╎    ╎    ╎    ╎    ╎   3      …olarData.jl:92; (::GRASS.var"#16#25"{…
     ╎    ╎    ╎    ╎    ╎    ╎    1      …tributes.jl:375; getindex(x::HDF5.At…
     ╎    ╎    ╎    ╎    ╎    ╎     1      …tributes.jl:66; open_attribute
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …tributes.jl:66; open_attribute
    1╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …unctions.jl:374; h5a_open(obj_id:…
     ╎    ╎    ╎    ╎    ╎    ╎    1      …eadwrite.jl:50; read(obj::HDF5.Attri…
     ╎    ╎    ╎    ╎    ╎    ╎     1      …versions.jl:238; get_jl_type
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …versions.jl:264; get_mem_compatibl…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …/helpers.jl:1004; h5t_get_native_…
    1╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …unctions.jl:7244; h5t_get_native…
     ╎    ╎    ╎    ╎    ╎    ╎    1      …eadwrite.jl:51; read(obj::HDF5.Attri…
     ╎    ╎    ╎    ╎    ╎    ╎     1      …eadwrite.jl:146; generic_read(::HDF…
    1╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …eadwrite.jl:178; _generic_read(::H…
     ╎    ╎    ╎    ╎    ╎    ╎   3      …olarData.jl:98; (::GRASS.var"#16#25"{…
     ╎    ╎    ╎    ╎    ╎    ╎    3      …se/array.jl:631; zeros
     ╎    ╎    ╎    ╎    ╎    ╎     3      …se/array.jl:633; zeros
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 3      …se/array.jl:637; zeros
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  3      …se/array.jl:395; fill!
    3╎    ╎    ╎    ╎    ╎    ╎    ╎   3      …se/array.jl:1021; setindex!
     ╎    ╎    ╎    ╎    ╎    ╎   1      …olarData.jl:99; (::GRASS.var"#16#25"{…
     ╎    ╎    ╎    ╎    ╎    ╎    1      …se/array.jl:631; zeros
     ╎    ╎    ╎    ╎    ╎    ╎     1      …se/array.jl:633; zeros
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …se/array.jl:637; zeros
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …se/array.jl:395; fill!
    1╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …se/array.jl:1021; setindex!
     ╎    ╎    ╎    ╎    ╎    ╎   1      …olarData.jl:100; (::GRASS.var"#16#25"…
     ╎    ╎    ╎    ╎    ╎    ╎    1      …se/array.jl:631; zeros
     ╎    ╎    ╎    ╎    ╎    ╎     1      …se/array.jl:633; zeros
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …se/array.jl:637; zeros
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …se/array.jl:395; fill!
    1╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …se/array.jl:1021; setindex!
    1╎    ╎    ╎    ╎    ╎    ╎   6      …olarData.jl:104; (::GRASS.var"#16#25"…
     ╎    ╎    ╎    ╎    ╎    ╎    1      …ighlevel.jl:74; getindex(parent::HDF…
     ╎    ╎    ╎    ╎    ╎    ╎     1      …ighlevel.jl:77; getindex(parent::HD…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …/objects.jl:37; open_object
    1╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …unctions.jl:2562; h5o_open(loc_id…
     ╎    ╎    ╎    ╎    ╎    ╎    2      …ighlevel.jl:74; getindex
     ╎    ╎    ╎    ╎    ╎    ╎     2      …ighlevel.jl:75; getindex(parent::HD…
    1╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …c/groups.jl:128; haskey
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …c/groups.jl:128; haskey(parent::HD…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …c/groups.jl:132; haskey(parent::H…
    1╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …unctions.jl:2076; h5l_exists(loc…
     ╎    ╎    ╎    ╎    ╎    ╎    1      …eadwrite.jl:50; read(obj::HDF5.Datas…
     ╎    ╎    ╎    ╎    ╎    ╎     1      …versions.jl:238; get_jl_type
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …versions.jl:274; get_mem_compatibl…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …unctions.jl:6974; h5t_close(dtype…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …ase/lock.jl:177; unlock
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …ase/lock.jl:182; (::Base.var"#_…
    1╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …ase/Base.jl:70; swapproperty!
     ╎    ╎    ╎    ╎    ╎    ╎    1      …eadwrite.jl:51; read(obj::HDF5.Datas…
     ╎    ╎    ╎    ╎    ╎    ╎     1      …eadwrite.jl:146; generic_read(::HDF…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …eadwrite.jl:180; _generic_read(::H…
    1╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …unctions.jl:796; h5d_read(dataset…
     ╎    ╎    ╎    ╎    ╎    ╎   5      …olarData.jl:105; (::GRASS.var"#16#25"…
     ╎    ╎    ╎    ╎    ╎    ╎    3      …ighlevel.jl:74; getindex(parent::HDF…
     ╎    ╎    ╎    ╎    ╎    ╎     3      …ighlevel.jl:77; getindex(parent::HD…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 3      …/objects.jl:37; open_object
    3╎    ╎    ╎    ╎    ╎    ╎    ╎  3      …unctions.jl:2562; h5o_open(loc_id…
     ╎    ╎    ╎    ╎    ╎    ╎    2      …eadwrite.jl:51; read(obj::HDF5.Datas…
     ╎    ╎    ╎    ╎    ╎    ╎     2      …eadwrite.jl:146; generic_read(::HDF…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 2      …eadwrite.jl:180; _generic_read(::H…
    2╎    ╎    ╎    ╎    ╎    ╎    ╎  2      …unctions.jl:796; h5d_read(dataset…
     ╎    ╎    ╎    ╎    ╎    ╎   3      …olarData.jl:106; (::GRASS.var"#16#25"…
     ╎    ╎    ╎    ╎    ╎    ╎    3      …eadwrite.jl:51; read(obj::HDF5.Datas…
     ╎    ╎    ╎    ╎    ╎    ╎     3      …eadwrite.jl:146; generic_read(::HDF…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 3      …eadwrite.jl:180; _generic_read(::H…
    3╎    ╎    ╎    ╎    ╎    ╎    ╎  3      …unctions.jl:796; h5d_read(dataset…
     ╎    ╎    ╎    ╎    ╎    ╎   4      …olarData.jl:107; (::GRASS.var"#16#25"…
     ╎    ╎    ╎    ╎    ╎    ╎    4      …eadwrite.jl:51; read(obj::HDF5.Datas…
     ╎    ╎    ╎    ╎    ╎    ╎     4      …eadwrite.jl:146; generic_read(::HDF…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …eadwrite.jl:164; _generic_read(::H…
    1╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …eadwrite.jl:387; _normalized_buff…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 3      …eadwrite.jl:180; _generic_read(::H…
    3╎    ╎    ╎    ╎    ╎    ╎    ╎  3      …unctions.jl:796; h5d_read(dataset…
     ╎    ╎    ╎    ╎    ╎    ╎   3      …olarData.jl:108; (::GRASS.var"#16#25"…
     ╎    ╎    ╎    ╎    ╎    ╎    1      …actarray.jl:1396; setindex!
     ╎    ╎    ╎    ╎    ╎    ╎     1      …ensional.jl:944; _setindex!
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …ensional.jl:955; _unsafe_setindex!…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …artesian.jl:64; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …ensional.jl:961; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …se/array.jl:945; iterate
    1╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …sentials.jl:13; getindex
     ╎    ╎    ╎    ╎    ╎    ╎    1      …se/array.jl:975; getindex
     ╎    ╎    ╎    ╎    ╎    ╎     1      …actarray.jl:831; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …se/array.jl:420; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …ase/boot.jl:486; Array
    1╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …ase/boot.jl:477; Array
     ╎    ╎    ╎    ╎    ╎    ╎    1      …eadwrite.jl:51; read(obj::HDF5.Datas…
     ╎    ╎    ╎    ╎    ╎    ╎     1      …eadwrite.jl:146; generic_read(::HDF…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …eadwrite.jl:180; _generic_read(::H…
    1╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …unctions.jl:796; h5d_read(dataset…
     ╎    ╎    ╎    ╎    ╎    ╎   12     …olarData.jl:114; (::GRASS.var"#16#25"…
     ╎    ╎    ╎    ╎    ╎    ╎    12     …/inputIO.jl:51; identify_bad_cols(bi…
     ╎    ╎    ╎    ╎    ╎    ╎     1      …/inputIO.jl:82; identify_bad_cols(b…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …atistics.jl:174; mean
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …atistics.jl:174; #mean#2
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …atistics.jl:187; _mean(f::typeof…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …educedim.jl:1011; sum
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …educedim.jl:1011; #sum#829
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …educedim.jl:1039; _sum
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …educedim.jl:1039; #_sum#855
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …educedim.jl:357; mapreduce
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …educedim.jl:357; #mapreduc…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …educedim.jl:371; _mapredu…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 1      …educedim.jl:324; mapreduc…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +2 1      …educedim.jl:316; _mapredu…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +3 1      …simdloop.jl:78; macro exp…
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +4 1      @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎    ╎     1      …/inputIO.jl:83; identify_bad_cols(b…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …atistics.jl:459; std
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …atistics.jl:459; #std#17
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …atistics.jl:467; _std(A::SubArra…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …atistics.jl:378; var
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …atistics.jl:378; #var#15
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …atistics.jl:382; _var
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …atistics.jl:174; mean
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …atistics.jl:174; #mean#2
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …atistics.jl:184; _mean(f::…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     1      @Base/set.jl:174; unique(i…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 1      @Base/set.jl:45; Set
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +2 1      …ase/dict.jl:70; Dict{Int6…
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +3 1      …ase/boot.jl:477; Array
     ╎    ╎    ╎    ╎    ╎    ╎     8      …/inputIO.jl:100; identify_bad_cols(…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 8      …rc/utils.jl:2; ismonotonic(A::SubA…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  3      …roadcast.jl:903; materialize
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   3      …roadcast.jl:928; copy
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …roadcast.jl:956; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …roadcast.jl:1014; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …roadcast.jl:1003; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …simdloop.jl:77; macro expans…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …roadcast.jl:1004; macro exp…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …bitarray.jl:702; setindex!
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …bitarray.jl:689; unsafe_b…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 1      …bitarray.jl:696; _unsafe_…
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +2 1      …se/array.jl:1021; setinde…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    2      …roadcast.jl:223; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     2      …roadcast.jl:226; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 2      …actarray.jl:876; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  2      …actarray.jl:877; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   2      …bitarray.jl:71; BitArray
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …bitarray.jl:37; BitArray
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …ase/boot.jl:477; Array
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …bitarray.jl:39; BitArray
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  5      …ensional.jl:979; diff
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   5      …ensional.jl:1011; diff
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    5      …ensional.jl:1019; diff(a::SubAr…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     4      …roadcast.jl:903; materialize
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 3      …roadcast.jl:928; copy
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  3      …roadcast.jl:223; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   3      …roadcast.jl:224; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    3      …actarray.jl:876; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     3      …actarray.jl:877; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 3      …ase/boot.jl:494; Array
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +2 3      …ase/boot.jl:486; Array
    3╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +3 3      …ase/boot.jl:477; Array
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …roadcast.jl:306; instantiate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …roadcast.jl:524; combine_axes
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …roadcast.jl:543; broadcast_…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …roadcast.jl:549; _bcs
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …roadcast.jl:555; _bcs1
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 1      …roadcast.jl:557; _bcsm
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …subarray.jl:184; view
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …actarray.jl:702; checkbounds
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …actarray.jl:687; checkbounds
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …actarray.jl:768; checkindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …actarray.jl:763; checkindex
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     1      @Base/int.jl:513; <
     ╎    ╎    ╎    ╎    ╎    ╎     1      …/inputIO.jl:105; identify_bad_cols(…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …roadcast.jl:903; materialize
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …roadcast.jl:928; copy
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …roadcast.jl:956; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …roadcast.jl:1014; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …roadcast.jl:1003; copyto!
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …simdloop.jl:0; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎     1      …/inputIO.jl:116; identify_bad_cols(…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …roadcast.jl:903; materialize
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …roadcast.jl:928; copy
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …roadcast.jl:956; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …roadcast.jl:1014; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …roadcast.jl:1000; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …roadcast.jl:983; preprocess
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …roadcast.jl:986; preprocess_…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …roadcast.jl:987; preprocess…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …roadcast.jl:983; preprocess
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …roadcast.jl:986; preproce…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 1      …roadcast.jl:987; preproce…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +2 1      …roadcast.jl:984; preproce…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +3 1      …roadcast.jl:977; broadcas…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +4 1      …actarray.jl:1481; unalias
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +5 1      …actarray.jl:1516; mightal…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +6 1      …actarray.jl:1539; dataids
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +7 1      …flection.jl:611; objectid
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +8 1      …flection.jl:617; _objectid
     ╎    ╎    ╎    ╎    ╎    ╎   1      …olarData.jl:116; (::GRASS.var"#16#25"…
     ╎    ╎    ╎    ╎    ╎    ╎    1      …rc/utils.jl:7; strip_columns(A::Vect…
     ╎    ╎    ╎    ╎    ╎    ╎     1      …actarray.jl:1291; getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …ensional.jl:889; _getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …ensional.jl:901; _unsafe_getindex…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …actarray.jl:831; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …se/array.jl:420; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …ase/boot.jl:486; Array
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …ase/boot.jl:477; Array
     ╎    ╎    ╎    ╎    ╎    ╎   1      …olarData.jl:124; (::GRASS.var"#16#25"…
     ╎    ╎    ╎    ╎    ╎    ╎    1      …educedim.jl:1010; minimum(a::Matrix{…
     ╎    ╎    ╎    ╎    ╎    ╎     1      …educedim.jl:1010; #minimum#840
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …educedim.jl:1014; _minimum
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …educedim.jl:1014; #_minimum#842
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …educedim.jl:1015; _minimum
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …educedim.jl:1015; #_minimum#843
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …educedim.jl:357; mapreduce
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …educedim.jl:357; #mapreduce#8…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …educedim.jl:365; _mapreduce_…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …e/reduce.jl:447; _mapreduce…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …e/reduce.jl:654; mapreduce…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …e/reduce.jl:632; _fast
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 1      …sentials.jl:647; ifelse
     ╎    ╎    ╎    ╎    ╎    ╎   161    …olarData.jl:125; (::GRASS.var"#16#25"…
     ╎    ╎    ╎    ╎    ╎    ╎    161    …/inputIO.jl:239; extrapolate_input_d…
    1╎    ╎    ╎    ╎    ╎    ╎     161    …/inputIO.jl:129; extrapolate_input_…
    1╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …/inputIO.jl:144; (::GRASS.var"#59#…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 13     …/inputIO.jl:138; extrapolate_input…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  13     …roadcast.jl:911; materialize!
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   13     …roadcast.jl:914; materialize!
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    13     …roadcast.jl:956; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     13     …roadcast.jl:1003; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 13     …simdloop.jl:77; macro expansi…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  13     …roadcast.jl:1004; macro expa…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   13     …roadcast.jl:636; getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    13     …roadcast.jl:682; _broadcas…
    3╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     13     …roadcast.jl:709; _broadca…
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 1      …erpolate.jl:0; (::GRASS.v…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 6      …erpolate.jl:11; (::GRASS.…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +2 6      …ase/sort.jl:292; searchso…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +3 6      …ase/sort.jl:292; #searchs…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +4 6      …ase/sort.jl:290; searchso…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +5 2      …ase/sort.jl:174; searchso…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +6 2      …perators.jl:276; !=
    2╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +7 2      …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +5 1      …ase/sort.jl:175; searchso…
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +6 1      @Base/int.jl:530; >>>
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +5 3      …ase/sort.jl:177; searchso…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +6 3      …ordering.jl:117; lt
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +7 3      …se/float.jl:549; isless
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +8 3      …se/float.jl:620; isnan
    3╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +9 3      …se/float.jl:535; !=
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 3      …erpolate.jl:14; (::GRASS.…
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +2 1      …se/float.jl:410; -
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +2 1      …se/float.jl:412; /
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +2 1      …subarray.jl:322; getindex
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +3 1      …actarray.jl:700; checkbou…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 15     …/inputIO.jl:144; extrapolate_input…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  15     …se/array.jl:2206; findfirst
    1╎    ╎    ╎    ╎    ╎    ╎    ╎   15     …se/array.jl:2155; findnext
   14╎    ╎    ╎    ╎    ╎    ╎    ╎    14     …/inputIO.jl:144; (::GRASS.var"#…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …/inputIO.jl:148; extrapolate_input…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …educedim.jl:1010; minimum
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …educedim.jl:1010; #minimum#840
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …educedim.jl:1014; _minimum
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …educedim.jl:1014; #_minimum#842
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …educedim.jl:1015; _minimum
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …educedim.jl:1015; #_minimum#…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …educedim.jl:357; mapreduce
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …educedim.jl:357; #mapreduc…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …educedim.jl:365; _mapredu…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 1      …e/reduce.jl:447; _mapredu…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +2 1      …e/reduce.jl:664; mapreduc…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +3 1      …subarray.jl:323; getindex
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +4 1      …sentials.jl:13; getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …/inputIO.jl:150; extrapolate_input…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …educedim.jl:1010; minimum
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …educedim.jl:1010; #minimum#840
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …educedim.jl:1014; _minimum
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …educedim.jl:1014; #_minimum#842
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …educedim.jl:1015; _minimum
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …educedim.jl:1015; #_minimum#…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …educedim.jl:357; mapreduce
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …educedim.jl:357; #mapreduc…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …educedim.jl:365; _mapredu…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 1      …e/reduce.jl:447; _mapredu…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +2 1      …e/reduce.jl:665; mapreduc…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +3 1      …ase/math.jl:899; min
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +4 1      …oatfuncs.jl:15; signbit
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 15     …/inputIO.jl:151; extrapolate_input…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  15     …se/array.jl:2206; findfirst
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   15     …se/array.jl:2155; findnext
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …subarray.jl:322; getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …actarray.jl:702; checkbounds
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …actarray.jl:687; checkbounds
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …actarray.jl:763; checkindex
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      @Base/int.jl:513; <
   13╎    ╎    ╎    ╎    ╎    ╎    ╎    14     …/inputIO.jl:151; (::GRASS.var"#…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …roadcast.jl:903; materialize(b…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …roadcast.jl:918; copy
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …roadcast.jl:636; getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …roadcast.jl:682; _broadcast…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …roadcast.jl:709; _broadcas…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …perators.jl:425; >=
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 1      …se/float.jl:537; <=
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 68     …/inputIO.jl:169; extrapolate_input…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  68     …c/common.jl:136; fit
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   68     …c/common.jl:136; #fit#21
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    68     …rd-basis.jl:554; fit
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     68     …rd-basis.jl:564; #fit#138
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 68     …c/common.jl:142; _fit
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  21     …c/common.jl:149; _fit(P::Typ…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   21     …rd-basis.jl:511; vander
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …rd-basis.jl:0; vander(P::T…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    2      …rd-basis.jl:516; vander(P:…
    2╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     2      …ase/boot.jl:479; Array
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    14     …rd-basis.jl:517; vander(P:…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     14     …se/array.jl:632; ones
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 14     …se/array.jl:636; ones
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +2 14     …ase/boot.jl:486; Array
   14╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +3 14     …ase/boot.jl:477; Array
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …rd-basis.jl:522; vander(P:…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …actarray.jl:1396; setinde…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 1      …ensional.jl:944; _setinde…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +2 1      …ensional.jl:955; _unsafe_…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +3 1      …artesian.jl:66; macro exp…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +4 1      …/indices.jl:393; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +5 1      …se/range.jl:901; iterate
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +6 1      …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    2      …rd-basis.jl:526; vander(P:…
    2╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     2      …se/array.jl:1021; setinde…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …rd-basis.jl:528; vander(P:…
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …se/range.jl:901; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  45     …c/common.jl:151; _fit(P::Typ…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   2      …c/common.jl:173; _wlstsq(va…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    2      …roadcast.jl:903; materiali…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     2      …roadcast.jl:928; copy
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 1      …roadcast.jl:956; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +2 1      …roadcast.jl:1003; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +3 1      …simdloop.jl:77; macro exp…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +4 1      …roadcast.jl:1004; macro e…
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +5 1      …se/array.jl:1021; setinde…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 1      …roadcast.jl:223; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +2 1      …roadcast.jl:224; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +3 1      …actarray.jl:876; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +4 1      …actarray.jl:877; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +5 1      …ase/boot.jl:494; Array
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +6 1      …ase/boot.jl:486; Array
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +7 1      …ase/boot.jl:477; Array
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   43     …c/common.jl:174; _wlstsq(va…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …diagonal.jl:285; *(D::Line…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …roadcast.jl:903; material…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 1      …roadcast.jl:928; copy
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +2 1      …roadcast.jl:223; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +3 1      …roadcast.jl:224; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +4 1      …actarray.jl:876; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +5 1      …actarray.jl:877; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +6 1      …ase/boot.jl:494; Array
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +7 1      …ase/boot.jl:486; Array
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +8 1      …ase/boot.jl:477; Array
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    2      …diagonal.jl:292; *(D::Line…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     2      …se/array.jl:418; similar
    2╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 2      …ase/boot.jl:479; Array
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    40     …/generic.jl:1126; \(A::Mat…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     30     …rAlgebra.jl:547; \
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 1      …rAlgebra.jl:569; ldiv(F::…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +2 1      …se/array.jl:388; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +3 1      …se/array.jl:368; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +4 1      …se/array.jl:374; _copyto_…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +5 1      …actarray.jl:702; checkbou…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +6 1      …actarray.jl:687; checkbou…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +7 1      …actarray.jl:768; checkind…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +8 1      …actarray.jl:763; checkind…
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +9 1      @Base/int.jl:86; -
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 28     …rAlgebra.jl:572; ldiv(F::…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +2 28     …a/src/qr.jl:610; ldiv!(A:…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +3 5      …rraymath.jl:41; vec
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +4 5      …pedarray.jl:117; reshape
    3╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +5 5      …pedarray.jl:51; reshape
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +6 2      …api/lock.jl:7; try_close_…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +7 2      …ase/lock.jl:238; trylock(…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +8 2      …api/lock.jl:8; #2
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +9 2      …/objects.jl:16; close
    2╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ +10 2      …unctions.jl:2176; h5o_clo…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +3 1      …pedarray.jl:117; reshape
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +4 1      …pedarray.jl:51; reshape
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +3 22     …a/src/qr.jl:612; ldiv!
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +4 1      …a/src/qr.jl:0; ldiv!(A::L…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +4 7      …a/src/qr.jl:564; ldiv!(A:…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +5 1      …c/lapack.jl:2519; laic1!(…
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +6 1      …ase/boot.jl:477; Array
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +5 6      …c/lapack.jl:2520; laic1!(…
    4╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +6 6      …ase/boot.jl:477; Array
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +7 2      …api/lock.jl:7; try_close_…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +8 2      …ase/lock.jl:238; trylock(…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +9 2      …api/lock.jl:8; #2
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ +10 2      …/objects.jl:16; close
    2╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ +11 2      …unctions.jl:2176; h5o_clo…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +4 3      …a/src/qr.jl:565; ldiv!(A:…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +5 1      …c/lapack.jl:2520; laic1!(…
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +6 1      …ase/boot.jl:477; Array
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +5 1      …c/lapack.jl:2521; laic1!(…
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +6 1      …ase/boot.jl:477; Array
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +5 1      …c/lapack.jl:2522; laic1!(…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +4 10     …a/src/qr.jl:589; ldiv!(A:…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +5 10     …bstractq.jl:347; lmul!
    3╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +6 3      …c/lapack.jl:2926; ormqr!(…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +6 7      …c/lapack.jl:2938; ormqr!(…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +7 7      …se/array.jl:1315; resize!
    7╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +8 7      …se/array.jl:1072; _growen…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +4 1      …a/src/qr.jl:590; ldiv!(A:…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +5 1      …iangular.jl:757; ldiv!
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +6 1      …iangular.jl:752; _ldiv!
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +7 1      …iangular.jl:830; generic_…
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +8 1      …c/lapack.jl:3557; trtrs!(…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 1      …rAlgebra.jl:577; ldiv(F::…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +2 1      …rAlgebra.jl:491; _cut_B
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +3 1      …se/array.jl:977; getindex
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +4 1      …se/array.jl:368; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     10     …a/src/qr.jl:417; qr
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 3      …a/src/qr.jl:419; #qr#222
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +2 3      …rAlgebra.jl:416; copy_sim…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +3 1      …se/array.jl:388; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +4 1      …se/array.jl:368; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +5 1      …se/array.jl:376; _copyto_…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +6 1      …se/array.jl:337; unsafe_c…
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +7 1      …ase/cmem.jl:26; memmove
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +3 2      …se/array.jl:420; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +4 2      …ase/boot.jl:487; Array
    2╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +5 2      …ase/boot.jl:479; Array
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 7      …a/src/qr.jl:420; #qr#222
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +2 7      …a/src/qr.jl:293; qr!
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +3 7      …c/lapack.jl:630; geqp3!(A…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +4 1      …se/array.jl:419; similar
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +5 1      …ase/boot.jl:477; Array
    6╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +4 6      …c/lapack.jl:426; geqp3!(A…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  2      …c/common.jl:157; _fit(P::Typ…
    2╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   2      …lynomial.jl:286; (Polynomia…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 32     …/inputIO.jl:170; extrapolate_input…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …actarray.jl:1291; getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …ensional.jl:889; _getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …ensional.jl:903; _unsafe_getind…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …ensional.jl:912; _unsafe_getin…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …artesian.jl:64; macro expansi…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …ensional.jl:917; macro expan…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …subarray.jl:323; getindex
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …sentials.jl:13; getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  31     …c/common.jl:136; fit
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   31     …c/common.jl:136; #fit#21
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    31     …rd-basis.jl:554; fit
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     31     …rd-basis.jl:564; #fit#138
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 31     …c/common.jl:142; _fit
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  29     …c/common.jl:153; _fit(P::Typ…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   29     …/generic.jl:1126; \(A::Matr…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    26     …rAlgebra.jl:547; \
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     26     …rAlgebra.jl:572; ldiv(F::…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 26     …a/src/qr.jl:610; ldiv!(A:…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +2 22     …rraymath.jl:41; vec
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +3 22     …pedarray.jl:117; reshape
   21╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +4 22     …pedarray.jl:51; reshape
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +5 1      …api/lock.jl:7; try_close_…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +6 1      …ase/lock.jl:238; trylock(…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +7 1      …api/lock.jl:8; #2
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +8 1      …/objects.jl:16; close
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +9 1      …unctions.jl:2176; h5o_clo…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +2 4      …a/src/qr.jl:612; ldiv!
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +3 1      …a/src/qr.jl:564; ldiv!(A:…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +4 1      …c/lapack.jl:2516; laic1!(…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +5 1      …actarray.jl:315; length
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +6 1      …subarray.jl:63; size
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +7 1      …subarray.jl:490; axes
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +8 1      …subarray.jl:495; _indices…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +9 1      …se/range.jl:706; axes
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ +10 1      …se/range.jl:469; oneto
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ +11 1      …se/range.jl:467; OneTo
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ +12 1      …se/range.jl:454; OneTo
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ +13 1      …romotion.jl:532; max
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ +14 1      …sentials.jl:647; ifelse
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +3 1      …a/src/qr.jl:565; ldiv!(A:…
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +4 1      …c/lapack.jl:2522; laic1!(…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +3 2      …a/src/qr.jl:589; ldiv!(A:…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +4 2      …bstractq.jl:347; lmul!
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +5 2      …c/lapack.jl:2938; ormqr!(…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +6 2      …se/array.jl:1315; resize!
    2╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +7 2      …se/array.jl:1072; _growen…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    3      …a/src/qr.jl:417; qr
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     3      …a/src/qr.jl:420; #qr#222
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 3      …a/src/qr.jl:293; qr!
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +2 3      …c/lapack.jl:630; geqp3!(A…
    3╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +3 3      …c/lapack.jl:426; geqp3!(A…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  2      …c/common.jl:157; _fit(P::Typ…
    2╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   2      …lynomial.jl:286; (Polynomia…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 7      …/inputIO.jl:173; extrapolate_input…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  5      …roadcast.jl:1341; broadcasted
    5╎    ╎    ╎    ╎    ╎    ╎    ╎   5      …roadcast.jl:1349; broadcasted
    2╎    ╎    ╎    ╎    ╎    ╎    ╎  2      …roadcast.jl:911; materialize!
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 6      …/inputIO.jl:174; extrapolate_input…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  6      …roadcast.jl:1341; broadcasted
    6╎    ╎    ╎    ╎    ╎    ╎    ╎   6      …roadcast.jl:1349; broadcasted
    1╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …erpolate.jl:2; (::GRASS.var"#f#7"{…
    1╎    ╎    ╎    ╎    ╎    ╎   1      …olarData.jl:149; (::GRASS.var"#16#25"…
     ╎    ╎    ╎    ╎    ╎    ╎   3      …olarData.jl:152; (::GRASS.var"#16#25"…
     ╎    ╎    ╎    ╎    ╎    ╎    3      …roadcast.jl:903; materialize(bc::Bas…
     ╎    ╎    ╎    ╎    ╎    ╎     3      …roadcast.jl:928; copy
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 3      …roadcast.jl:956; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  3      …roadcast.jl:1003; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   3      …simdloop.jl:77; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    3      …roadcast.jl:1004; macro expansi…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     2      …roadcast.jl:636; getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 2      …roadcast.jl:682; _broadcast_g…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  2      …roadcast.jl:709; _broadcast_…
    2╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   2      …se/float.jl:409; +
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …ensional.jl:698; setindex!
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …se/array.jl:1024; setindex!
     ╎    ╎    ╎    ╎    ╎    ╎  2      …src/file.jl:103; (::HDF5.var"#17#18"{H…
     ╎    ╎    ╎    ╎    ╎    ╎   2      …src/file.jl:165; close
    2╎    ╎    ╎    ╎    ╎    ╎    2      …unctions.jl:1068; h5f_close(file_id:…
     ╎    ╎    ╎    ╎    ╎   1      …olarData.jl:169; SolarData(fname::String; …
    1╎    ╎    ╎    ╎    ╎    1      …subarray.jl:181; view(::Matrix{Float64}, …
     ╎    ╎    ╎    ╎    132084 …ce_eclipse.jl:62; synth_Eclipse_cpu(spec::Spec…
     ╎    ╎    ╎    ╎     132084 …im_eclipse.jl:1; kwcall(::@NamedTuple{skip_ti…
   99╎    ╎    ╎    ╎    ╎ 132084 …m_eclipse.jl:11; disk_sim_eclipse(spec::Spec…
    1╎    ╎    ╎    ╎    ╎  1      …se/array.jl:851; _collect(c::Matrix{Float64…
   31╎    ╎    ╎    ╎    ╎  31     …se/array.jl:371; _copyto_impl!(dest::Vector…
   15╎    ╎    ╎    ╎    ╎  15     …se/array.jl:2031; vcat(::Vector{Float64}, :…
    1╎    ╎    ╎    ╎    ╎  1      …rraymath.jl:24; /(A::Matrix{Float64}, B::Fl…
    2╎    ╎    ╎    ╎    ╎  2      …bitarray.jl:28; BitMatrix(::UndefInitialize…
    1╎    ╎    ╎    ╎    ╎  1      …bitarray.jl:283; copy_to_bitarray_chunks!(B…
    5╎    ╎    ╎    ╎    ╎  5      …/generic.jl:0; generic_norm2(x::SubArray{Fl…
    3╎    ╎    ╎    ╎    ╎  3      …/generic.jl:462; generic_norm2(x::SubArray{…
    7╎    ╎    ╎    ╎    ╎  7      …/generic.jl:595; norm(itr::Vector{Float64},…
    1╎    ╎    ╎    ╎    ╎  1      …c/matmul.jl:0; gemv!(y::Vector{Float64}, tA…
   26╎    ╎    ╎    ╎    ╎  26     …c/matmul.jl:401; gemv!(y::Vector{Float64}, …
   80╎    ╎    ╎    ╎    ╎  80     …pse_comp.jl:0; eclipse_compute_quantities!(…
   11╎    ╎    ╎    ╎    ╎  11     …pse_comp.jl:46; eclipse_compute_quantities!…
   18╎    ╎    ╎    ╎    ╎  22     …pse_comp.jl:48; eclipse_compute_quantities!…
    2╎    ╎    ╎    ╎    ╎   2      …se/array.jl:1024; setindex!(::Matrix{Float…
    1╎    ╎    ╎    ╎    ╎   1      …se/range.jl:952; getindex(r::StepRangeLen{…
     ╎    ╎    ╎    ╎    ╎   1      …se/range.jl:956; getindex(r::StepRangeLen{…
     ╎    ╎    ╎    ╎    ╎    1      …recision.jl:483; unsafe_getindex
     ╎    ╎    ╎    ╎    ╎     1      …recision.jl:84; add12
    1╎    ╎    ╎    ╎    ╎    ╎ 1      …se/float.jl:610; abs
    5╎    ╎    ╎    ╎    ╎  14     …pse_comp.jl:49; eclipse_compute_quantities!…
    4╎    ╎    ╎    ╎    ╎   4      …se/array.jl:1024; setindex!(::Matrix{Float…
    5╎    ╎    ╎    ╎    ╎   5      …sentials.jl:14; getindex(::Matrix{Float64}…
   32╎    ╎    ╎    ╎    ╎  59     …pse_comp.jl:52; eclipse_compute_quantities!…
    1╎    ╎    ╎    ╎    ╎   1      @Base/int.jl:87; +(x::Int64, y::Int64)
     ╎    ╎    ╎    ╎    ╎   26     …se/range.jl:149; kwcall(::@NamedTuple{leng…
     ╎    ╎    ╎    ╎    ╎    26     …se/range.jl:149; #range#76
    1╎    ╎    ╎    ╎    ╎     26     …se/range.jl:166; _range
    1╎    ╎    ╎    ╎    ╎    ╎ 1      …recision.jl:669; _linspace(start::Float…
    2╎    ╎    ╎    ╎    ╎    ╎ 2      …recision.jl:774; lcm_unchecked(a::Int64…
    1╎    ╎    ╎    ╎    ╎    ╎ 1      …recision.jl:646; range_start_stop_lengt…
     ╎    ╎    ╎    ╎    ╎    ╎ 9      …recision.jl:653; range_start_stop_lengt…
    1╎    ╎    ╎    ╎    ╎    ╎  1      …recision.jl:756; rat(x::Float64)
     ╎    ╎    ╎    ╎    ╎    ╎  3      …recision.jl:761; rat(x::Float64)
    1╎    ╎    ╎    ╎    ╎    ╎   1      …se/float.jl:610; abs
     ╎    ╎    ╎    ╎    ╎    ╎   2      …romotion.jl:463; <=
    2╎    ╎    ╎    ╎    ╎    ╎    2      …se/float.jl:537; <=
     ╎    ╎    ╎    ╎    ╎    ╎  1      …recision.jl:762; rat(x::Float64)
    1╎    ╎    ╎    ╎    ╎    ╎   1      …se/float.jl:902; trunc
     ╎    ╎    ╎    ╎    ╎    ╎  1      …recision.jl:764; rat(x::Float64)
    1╎    ╎    ╎    ╎    ╎    ╎   1      @Base/int.jl:88; *
     ╎    ╎    ╎    ╎    ╎    ╎  3      …recision.jl:767; rat(x::Float64)
    3╎    ╎    ╎    ╎    ╎    ╎   3      …se/float.jl:534; ==
     ╎    ╎    ╎    ╎    ╎    ╎ 3      …recision.jl:654; range_start_stop_lengt…
    1╎    ╎    ╎    ╎    ╎    ╎  1      …recision.jl:756; rat(x::Float64)
     ╎    ╎    ╎    ╎    ╎    ╎  1      …recision.jl:762; rat(x::Float64)
    1╎    ╎    ╎    ╎    ╎    ╎   1      …se/float.jl:902; trunc
     ╎    ╎    ╎    ╎    ╎    ╎  1      …recision.jl:767; rat(x::Float64)
    1╎    ╎    ╎    ╎    ╎    ╎   1      …se/float.jl:412; /
     ╎    ╎    ╎    ╎    ╎    ╎ 2      …recision.jl:656; range_start_stop_lengt…
     ╎    ╎    ╎    ╎    ╎    ╎  2      …recision.jl:774; lcm_unchecked(a::Int6…
    1╎    ╎    ╎    ╎    ╎    ╎   1      @Base/int.jl:88; *
     ╎    ╎    ╎    ╎    ╎    ╎   1      …intfuncs.jl:60; gcd
     ╎    ╎    ╎    ╎    ╎    ╎    1      …intfuncs.jl:84; _gcd
     ╎    ╎    ╎    ╎    ╎    ╎     1      @Base/int.jl:534; >>
    1╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      @Base/int.jl:527; >>
     ╎    ╎    ╎    ╎    ╎    ╎ 3      …recision.jl:666; range_start_stop_lengt…
    1╎    ╎    ╎    ╎    ╎    ╎  1      …recision.jl:669; _linspace(start::Floa…
     ╎    ╎    ╎    ╎    ╎    ╎  1      …recision.jl:680; _linspace(start::Floa…
     ╎    ╎    ╎    ╎    ╎    ╎   1      …se/float.jl:385; round
     ╎    ╎    ╎    ╎    ╎    ╎    1      …se/float.jl:902; trunc
    1╎    ╎    ╎    ╎    ╎    ╎     1      …se/float.jl:537; <=
     ╎    ╎    ╎    ╎    ╎    ╎  1      …recision.jl:703; _linspace(start::Floa…
    1╎    ╎    ╎    ╎    ╎    ╎   1      …se/float.jl:409; +
    4╎    ╎    ╎    ╎    ╎    ╎ 4      …recision.jl:756; rat(x::Float64)
   23╎    ╎    ╎    ╎    ╎  35     …pse_comp.jl:53; eclipse_compute_quantities!…
    5╎    ╎    ╎    ╎    ╎   5      …sentials.jl:14; getindex(::Matrix{Float64}…
     ╎    ╎    ╎    ╎    ╎   7      …se/range.jl:149; kwcall(::@NamedTuple{leng…
     ╎    ╎    ╎    ╎    ╎    7      …se/range.jl:149; #range#76
     ╎    ╎    ╎    ╎    ╎     7      …se/range.jl:166; _range
     ╎    ╎    ╎    ╎    ╎    ╎ 4      …recision.jl:653; range_start_stop_lengt…
    1╎    ╎    ╎    ╎    ╎    ╎  1      …recision.jl:0; rat(x::Float64)
     ╎    ╎    ╎    ╎    ╎    ╎  1      …recision.jl:766; rat(x::Float64)
     ╎    ╎    ╎    ╎    ╎    ╎   1      @Base/int.jl:188; abs
    1╎    ╎    ╎    ╎    ╎    ╎    1      @Base/int.jl:142; flipsign
    1╎    ╎    ╎    ╎    ╎    ╎  2      …recision.jl:767; rat(x::Float64)
    1╎    ╎    ╎    ╎    ╎    ╎   1      …se/float.jl:534; ==
     ╎    ╎    ╎    ╎    ╎    ╎ 2      …recision.jl:654; range_start_stop_lengt…
     ╎    ╎    ╎    ╎    ╎    ╎  2      …recision.jl:767; rat(x::Float64)
    2╎    ╎    ╎    ╎    ╎    ╎   2      …se/float.jl:534; ==
     ╎    ╎    ╎    ╎    ╎    ╎ 1      …recision.jl:656; range_start_stop_lengt…
     ╎    ╎    ╎    ╎    ╎    ╎  1      …recision.jl:774; lcm_unchecked(a::Int6…
    1╎    ╎    ╎    ╎    ╎    ╎   1      @Base/int.jl:88; *
    9╎    ╎    ╎    ╎    ╎  26     …pse_comp.jl:54; eclipse_compute_quantities!…
     ╎    ╎    ╎    ╎    ╎   1      …geometry.jl:12; get_grid_centers(grid::Ste…
     ╎    ╎    ╎    ╎    ╎    1      …se/range.jl:836; first
     ╎    ╎    ╎    ╎    ╎     1      …recision.jl:483; unsafe_getindex
     ╎    ╎    ╎    ╎    ╎    ╎ 1      …recision.jl:84; add12
    1╎    ╎    ╎    ╎    ╎    ╎  1      …se/float.jl:610; abs
     ╎    ╎    ╎    ╎    ╎   16     …geometry.jl:14; get_grid_centers(grid::Ste…
     ╎    ╎    ╎    ╎    ╎    16     …se/range.jl:149; range
     ╎    ╎    ╎    ╎    ╎     16     …se/range.jl:149; #range#76
     ╎    ╎    ╎    ╎    ╎    ╎ 16     …se/range.jl:166; _range
     ╎    ╎    ╎    ╎    ╎    ╎  4      …recision.jl:653; range_start_stop_leng…
    1╎    ╎    ╎    ╎    ╎    ╎   4      …recision.jl:767; rat(x::Float64)
     ╎    ╎    ╎    ╎    ╎    ╎    2      …sentials.jl:522; oftype
     ╎    ╎    ╎    ╎    ╎    ╎     2      …e/number.jl:7; convert
    2╎    ╎    ╎    ╎    ╎    ╎    ╎ 2      …se/float.jl:159; Float64
    1╎    ╎    ╎    ╎    ╎    ╎    1      …se/float.jl:534; ==
     ╎    ╎    ╎    ╎    ╎    ╎  8      …recision.jl:654; range_start_stop_leng…
     ╎    ╎    ╎    ╎    ╎    ╎   1      …recision.jl:762; rat(x::Float64)
    1╎    ╎    ╎    ╎    ╎    ╎    1      …se/float.jl:902; trunc
     ╎    ╎    ╎    ╎    ╎    ╎   1      …recision.jl:764; rat(x::Float64)
    1╎    ╎    ╎    ╎    ╎    ╎    1      @Base/int.jl:88; *
     ╎    ╎    ╎    ╎    ╎    ╎   3      …recision.jl:766; rat(x::Float64)
     ╎    ╎    ╎    ╎    ╎    ╎    2      @Base/int.jl:188; abs
    2╎    ╎    ╎    ╎    ╎    ╎     2      @Base/int.jl:142; flipsign
     ╎    ╎    ╎    ╎    ╎    ╎    1      …romotion.jl:532; max
    1╎    ╎    ╎    ╎    ╎    ╎     1      …sentials.jl:647; ifelse
    1╎    ╎    ╎    ╎    ╎    ╎   3      …recision.jl:767; rat(x::Float64)
    2╎    ╎    ╎    ╎    ╎    ╎    2      …se/float.jl:534; ==
     ╎    ╎    ╎    ╎    ╎    ╎  2      …recision.jl:656; range_start_stop_leng…
     ╎    ╎    ╎    ╎    ╎    ╎   2      …recision.jl:774; lcm_unchecked(a::Int…
    1╎    ╎    ╎    ╎    ╎    ╎    1      @Base/int.jl:88; *
    1╎    ╎    ╎    ╎    ╎    ╎    1      @Base/int.jl:295; div
     ╎    ╎    ╎    ╎    ╎    ╎  1      …recision.jl:661; range_start_stop_leng…
     ╎    ╎    ╎    ╎    ╎    ╎   1      @Base/int.jl:97; /
    1╎    ╎    ╎    ╎    ╎    ╎    1      …se/float.jl:412; /
     ╎    ╎    ╎    ╎    ╎    ╎  1      …recision.jl:666; range_start_stop_leng…
     ╎    ╎    ╎    ╎    ╎    ╎   1      …recision.jl:703; _linspace(start::Flo…
     ╎    ╎    ╎    ╎    ╎    ╎    1      …ase/math.jl:908; max
    1╎    ╎    ╎    ╎    ╎    ╎     1      …se/float.jl:410; -
    2╎    ╎    ╎    ╎    ╎  17     …pse_comp.jl:55; eclipse_compute_quantities!…
     ╎    ╎    ╎    ╎    ╎   15     …geometry.jl:14; get_grid_centers(grid::Ste…
     ╎    ╎    ╎    ╎    ╎    15     …se/range.jl:149; range
     ╎    ╎    ╎    ╎    ╎     15     …se/range.jl:149; #range#76
     ╎    ╎    ╎    ╎    ╎    ╎ 15     …se/range.jl:166; _range
     ╎    ╎    ╎    ╎    ╎    ╎  6      …recision.jl:653; range_start_stop_leng…
     ╎    ╎    ╎    ╎    ╎    ╎   2      …recision.jl:761; rat(x::Float64)
     ╎    ╎    ╎    ╎    ╎    ╎    2      …romotion.jl:463; <=
    2╎    ╎    ╎    ╎    ╎    ╎     2      …se/float.jl:537; <=
     ╎    ╎    ╎    ╎    ╎    ╎   1      …recision.jl:762; rat(x::Float64)
    1╎    ╎    ╎    ╎    ╎    ╎    1      …se/float.jl:902; trunc
     ╎    ╎    ╎    ╎    ╎    ╎   2      …recision.jl:766; rat(x::Float64)
     ╎    ╎    ╎    ╎    ╎    ╎    1      @Base/int.jl:188; abs
    1╎    ╎    ╎    ╎    ╎    ╎     1      @Base/int.jl:142; flipsign
     ╎    ╎    ╎    ╎    ╎    ╎    1      …romotion.jl:532; max
    1╎    ╎    ╎    ╎    ╎    ╎     1      …sentials.jl:647; ifelse
     ╎    ╎    ╎    ╎    ╎    ╎   1      …recision.jl:767; rat(x::Float64)
     ╎    ╎    ╎    ╎    ╎    ╎    1      …sentials.jl:522; oftype
     ╎    ╎    ╎    ╎    ╎    ╎     1      …e/number.jl:7; convert
    1╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …se/float.jl:159; Float64
     ╎    ╎    ╎    ╎    ╎    ╎  8      …recision.jl:654; range_start_stop_leng…
    1╎    ╎    ╎    ╎    ╎    ╎   1      …recision.jl:0; rat(x::Float64)
    1╎    ╎    ╎    ╎    ╎    ╎   1      …recision.jl:756; rat(x::Float64)
     ╎    ╎    ╎    ╎    ╎    ╎   1      …recision.jl:761; rat(x::Float64)
     ╎    ╎    ╎    ╎    ╎    ╎    1      …romotion.jl:463; <=
    1╎    ╎    ╎    ╎    ╎    ╎     1      …se/float.jl:537; <=
     ╎    ╎    ╎    ╎    ╎    ╎   1      …recision.jl:764; rat(x::Float64)
    1╎    ╎    ╎    ╎    ╎    ╎    1      @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎    ╎   2      …recision.jl:766; rat(x::Float64)
    2╎    ╎    ╎    ╎    ╎    ╎    2      @Base/int.jl:514; <=
     ╎    ╎    ╎    ╎    ╎    ╎   2      …recision.jl:767; rat(x::Float64)
    1╎    ╎    ╎    ╎    ╎    ╎    1      …se/float.jl:412; /
    1╎    ╎    ╎    ╎    ╎    ╎    1      …se/float.jl:534; ==
     ╎    ╎    ╎    ╎    ╎    ╎  1      …recision.jl:666; range_start_stop_leng…
     ╎    ╎    ╎    ╎    ╎    ╎   1      …recision.jl:711; _linspace(start::Flo…
     ╎    ╎    ╎    ╎    ╎    ╎    1      …recision.jl:362; steprangelen_hp
     ╎    ╎    ╎    ╎    ╎    ╎     1      …recision.jl:375; StepRangeLen
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …se/range.jl:508; StepRangeLen
    1╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …se/range.jl:0; StepRangeLen
   20╎    ╎    ╎    ╎    ╎  72     …pse_comp.jl:56; eclipse_compute_quantities!…
   52╎    ╎    ╎    ╎    ╎   52     …terators.jl:1024; product(::StepRangeLen{F…
   64╎    ╎    ╎    ╎    ╎  59900  …pse_comp.jl:59; eclipse_compute_quantities!…
   32╎    ╎    ╎    ╎    ╎   59803  …actarray.jl:3313; map(f::Function, A::Base…
    1╎    ╎    ╎    ╎    ╎    1      …se/array.jl:827; collect(itr::Base.Genera…
     ╎    ╎    ╎    ╎    ╎    166    …se/array.jl:834; collect(itr::Base.Genera…
     ╎    ╎    ╎    ╎    ╎     2      …enerator.jl:44; iterate
     ╎    ╎    ╎    ╎    ╎    ╎ 2      …terators.jl:1102; iterate
     ╎    ╎    ╎    ╎    ╎    ╎  2      …terators.jl:1094; _piterate
     ╎    ╎    ╎    ╎    ╎    ╎   2      …se/range.jl:892; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    2      …se/range.jl:894; iterate
     ╎    ╎    ╎    ╎    ╎    ╎     2      …recision.jl:481; unsafe_getindex
    2╎    ╎    ╎    ╎    ╎    ╎    ╎ 2      @Base/int.jl:86; -
     ╎    ╎    ╎    ╎    ╎     164    …enerator.jl:47; iterate
    9╎    ╎    ╎    ╎    ╎    ╎ 164    …pse_comp.jl:59; (::GRASS.var"#283#287")…
     ╎    ╎    ╎    ╎    ╎    ╎  5      …CE/src/p.jl:282; pgrrec(body::String, …
     ╎    ╎    ╎    ╎    ╎    ╎   5      …ase/boot.jl:491; Array
    5╎    ╎    ╎    ╎    ╎    ╎    5      …ase/boot.jl:477; Array
  146╎    ╎    ╎    ╎    ╎    ╎  147    …CE/src/p.jl:283; pgrrec(body::String, …
     ╎    ╎    ╎    ╎    ╎    ╎   1      @Base/c.jl:216; unsafe_convert
    1╎    ╎    ╎    ╎    ╎    ╎    1      @Base/c.jl:209; containsnul
    1╎    ╎    ╎    ╎    ╎    ╎  3      …CE/src/p.jl:286; pgrrec(body::String, …
    2╎    ╎    ╎    ╎    ╎    ╎   2      …rc/SPICE.jl:22; handleerror()
     ╎    ╎    ╎    ╎    ╎    57     …se/array.jl:839; collect(itr::Base.Genera…
     ╎    ╎    ╎    ╎    ╎     57     …se/array.jl:723; _array_for
     ╎    ╎    ╎    ╎    ╎    ╎ 57     …actarray.jl:876; similar
     ╎    ╎    ╎    ╎    ╎    ╎  57     …actarray.jl:877; similar
     ╎    ╎    ╎    ╎    ╎    ╎   57     …ase/boot.jl:487; Array
   57╎    ╎    ╎    ╎    ╎    ╎    57     …ase/boot.jl:479; Array
     ╎    ╎    ╎    ╎    ╎    59545  …se/array.jl:844; collect(itr::Base.Genera…
     ╎    ╎    ╎    ╎    ╎     59545  …se/array.jl:870; collect_to_with_first!
    5╎    ╎    ╎    ╎    ╎    ╎ 5      …se/array.jl:0; collect_to!(dest::Matrix…
    4╎    ╎    ╎    ╎    ╎    ╎ 4      …se/array.jl:887; collect_to!(dest::Matr…
     ╎    ╎    ╎    ╎    ╎    ╎ 59432  …se/array.jl:892; collect_to!(dest::Matr…
   19╎    ╎    ╎    ╎    ╎    ╎  484    …enerator.jl:44; iterate
     ╎    ╎    ╎    ╎    ╎    ╎   465    …terators.jl:1123; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    440    …terators.jl:1110; _piterate1
     ╎    ╎    ╎    ╎    ╎    ╎     4      …se/range.jl:892; iterate
    4╎    ╎    ╎    ╎    ╎    ╎    ╎ 4      @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎    ╎     436    …se/range.jl:894; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 8      …recision.jl:481; unsafe_getindex
    8╎    ╎    ╎    ╎    ╎    ╎    ╎  8      @Base/int.jl:86; -
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 60     …recision.jl:482; unsafe_getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  60     …romotion.jl:423; *
   52╎    ╎    ╎    ╎    ╎    ╎    ╎   52     …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   8      …romotion.jl:393; promote
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    8      …romotion.jl:370; _promote
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     8      …e/number.jl:7; convert
    8╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 8      …se/float.jl:159; Float64
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 210    …recision.jl:483; unsafe_getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  120    …recision.jl:84; add12
   58╎    ╎    ╎    ╎    ╎    ╎    ╎   58     …sentials.jl:647; ifelse
    5╎    ╎    ╎    ╎    ╎    ╎    ╎   5      …se/float.jl:610; abs
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   57     …perators.jl:378; >
   57╎    ╎    ╎    ╎    ╎    ╎    ╎    57     …se/float.jl:536; <
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  90     …recision.jl:85; add12
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   28     …recision.jl:49; canonicalize2
   28╎    ╎    ╎    ╎    ╎    ╎    ╎    28     …se/float.jl:409; +
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   62     …recision.jl:50; canonicalize2
   49╎    ╎    ╎    ╎    ╎    ╎    ╎    49     …se/float.jl:409; +
   13╎    ╎    ╎    ╎    ╎    ╎    ╎    13     …se/float.jl:410; -
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 158    …recision.jl:484; unsafe_getindex
  158╎    ╎    ╎    ╎    ╎    ╎    ╎  158    …se/float.jl:409; +
     ╎    ╎    ╎    ╎    ╎    ╎    25     …terators.jl:1114; _piterate1
     ╎    ╎    ╎    ╎    ╎    ╎     25     …terators.jl:1110; _piterate1
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 2      …se/range.jl:892; iterate
    2╎    ╎    ╎    ╎    ╎    ╎    ╎  2      @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 23     …se/range.jl:894; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  8      …recision.jl:482; unsafe_getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   8      …romotion.jl:423; *
    8╎    ╎    ╎    ╎    ╎    ╎    ╎    8      …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  11     …recision.jl:483; unsafe_getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   9      …recision.jl:84; add12
    2╎    ╎    ╎    ╎    ╎    ╎    ╎    2      …sentials.jl:647; ifelse
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    7      …perators.jl:378; >
    7╎    ╎    ╎    ╎    ╎    ╎    ╎     7      …se/float.jl:536; <
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   2      …recision.jl:85; add12
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    2      …recision.jl:50; canonicalize2
    2╎    ╎    ╎    ╎    ╎    ╎    ╎     2      …se/float.jl:409; +
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  4      …recision.jl:484; unsafe_getindex
    4╎    ╎    ╎    ╎    ╎    ╎    ╎   4      …se/float.jl:409; +
   71╎    ╎    ╎    ╎    ╎    ╎  58948  …enerator.jl:47; iterate
 1712╎    ╎    ╎    ╎    ╎    ╎   58877  …pse_comp.jl:59; (::GRASS.var"#283#287…
   91╎    ╎    ╎    ╎    ╎    ╎    91     …rc/SPICE.jl:21; handleerror()
   17╎    ╎    ╎    ╎    ╎    ╎    17     …rc/SPICE.jl:30; handleerror()
 1670╎    ╎    ╎    ╎    ╎    ╎    1670   …CE/src/p.jl:281; pgrrec(body::String…
     ╎    ╎    ╎    ╎    ╎    ╎    2548   …CE/src/p.jl:282; pgrrec(body::String…
     ╎    ╎    ╎    ╎    ╎    ╎     2548   …ase/boot.jl:491; Array
 2548╎    ╎    ╎    ╎    ╎    ╎    ╎ 2548   …ase/boot.jl:477; Array
52232╎    ╎    ╎    ╎    ╎    ╎    52330  …CE/src/p.jl:283; pgrrec(body::String…
     ╎    ╎    ╎    ╎    ╎    ╎     98     @Base/c.jl:216; unsafe_convert
   13╎    ╎    ╎    ╎    ╎    ╎    ╎ 13     …ase/Base.jl:244; sizeof
   84╎    ╎    ╎    ╎    ╎    ╎    ╎ 85     @Base/c.jl:209; containsnul
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …perators.jl:276; !=
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …/pointer.jl:278; ==
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …romotion.jl:521; ==
  120╎    ╎    ╎    ╎    ╎    ╎    509    …CE/src/p.jl:286; pgrrec(body::String…
   14╎    ╎    ╎    ╎    ╎    ╎     14     …rc/SPICE.jl:0; handleerror()
  180╎    ╎    ╎    ╎    ╎    ╎     180    …rc/SPICE.jl:21; handleerror()
  180╎    ╎    ╎    ╎    ╎    ╎     180    …rc/SPICE.jl:22; handleerror()
   15╎    ╎    ╎    ╎    ╎    ╎     15     …rc/SPICE.jl:30; handleerror()
     ╎    ╎    ╎    ╎    ╎    ╎ 94     …se/array.jl:896; collect_to!(dest::Matr…
   94╎    ╎    ╎    ╎    ╎    ╎  94     …se/array.jl:1021; setindex!
    3╎    ╎    ╎    ╎    ╎    ╎ 3      …se/array.jl:902; collect_to!(dest::Matr…
    7╎    ╎    ╎    ╎    ╎    ╎ 7      …pse_comp.jl:59; (::GRASS.var"#283#287")…
    1╎    ╎    ╎    ╎    ╎    1      …ase/boot.jl:0; collect(itr::Base.Generato…
    1╎    ╎    ╎    ╎    ╎    1      …enerator.jl:32; Base.Generator(f::GRASS.v…
    1╎    ╎    ╎    ╎    ╎   1      …se/array.jl:371; _copyto_impl!(dest::Matri…
    2╎    ╎    ╎    ╎    ╎   2      …roadcast.jl:910; materialize!(dest::Matrix…
     ╎    ╎    ╎    ╎    ╎   30     …roadcast.jl:911; materialize!(dest::Matrix…
     ╎    ╎    ╎    ╎    ╎    30     …roadcast.jl:914; materialize!
     ╎    ╎    ╎    ╎    ╎     1      …actarray.jl:98; axes
    1╎    ╎    ╎    ╎    ╎    ╎ 1      …se/array.jl:192; size
     ╎    ╎    ╎    ╎    ╎     29     …roadcast.jl:956; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎ 29     …roadcast.jl:997; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎  29     …se/array.jl:388; copyto!
    7╎    ╎    ╎    ╎    ╎    ╎   29     …se/array.jl:368; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎    22     …se/array.jl:376; _copyto_impl!(dest:…
   22╎    ╎    ╎    ╎    ╎    ╎     22     …se/array.jl:334; unsafe_copyto!
   43╎    ╎    ╎    ╎    ╎  5838   …pse_comp.jl:62; eclipse_compute_quantities!…
   53╎    ╎    ╎    ╎    ╎   5761   …actarray.jl:3313; map(f::Function, A::Base…
     ╎    ╎    ╎    ╎    ╎    39     …se/array.jl:834; collect(itr::Base.Genera…
     ╎    ╎    ╎    ╎    ╎     39     …enerator.jl:47; iterate
     ╎    ╎    ╎    ╎    ╎    ╎ 39     …pse_comp.jl:62; #284
    3╎    ╎    ╎    ╎    ╎    ╎  3      …ial/trig.jl:98; cos(x::Float64)
    9╎    ╎    ╎    ╎    ╎    ╎  36     …_physics.jl:57; v_scalar(lat::Float64,…
    1╎    ╎    ╎    ╎    ╎    ╎   1      …se/float.jl:412; /(x::Float64, y::Flo…
   16╎    ╎    ╎    ╎    ╎    ╎   16     …perators.jl:587; *
    2╎    ╎    ╎    ╎    ╎    ╎   2      …ial/trig.jl:0; cos(x::Float64)
     ╎    ╎    ╎    ╎    ╎    ╎   1      …ial/trig.jl:99; cos(x::Float64)
    1╎    ╎    ╎    ╎    ╎    ╎    1      …se/float.jl:610; abs
     ╎    ╎    ╎    ╎    ╎    ╎   1      …ial/trig.jl:101; cos(x::Float64)
    1╎    ╎    ╎    ╎    ╎    ╎    1      …se/float.jl:536; <
     ╎    ╎    ╎    ╎    ╎    ╎   6      …_physics.jl:50; rotation_period
     ╎    ╎    ╎    ╎    ╎    ╎    1      …_physics.jl:52; rotation_period(ϕ::F…
     ╎    ╎    ╎    ╎    ╎    ╎     1      …ial/trig.jl:35; sin(x::Float64)
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …ial/trig.jl:81; sin_kernel
    1╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …se/float.jl:409; +
     ╎    ╎    ╎    ╎    ╎    ╎    5      …_physics.jl:53; rotation_period(ϕ::F…
    1╎    ╎    ╎    ╎    ╎    ╎     1      …se/float.jl:412; /
    1╎    ╎    ╎    ╎    ╎    ╎     4      …ase/math.jl:1204; ^
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 2      …ase/math.jl:1248; ^(x::Float64, n:…
    1╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …ase/math.jl:0; pow_body(x::Float6…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …ase/math.jl:1263; pow_body(x::Flo…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …perators.jl:378; >
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    1      @Base/int.jl:83; <
    1╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …ase/math.jl:1251; pow_body(x::Floa…
   38╎    ╎    ╎    ╎    ╎    73     …se/array.jl:839; collect(itr::Base.Genera…
     ╎    ╎    ╎    ╎    ╎     35     …se/array.jl:723; _array_for(::Type{Float…
     ╎    ╎    ╎    ╎    ╎    ╎ 35     …actarray.jl:876; similar
     ╎    ╎    ╎    ╎    ╎    ╎  35     …actarray.jl:877; similar
     ╎    ╎    ╎    ╎    ╎    ╎   35     …ase/boot.jl:487; Array
   35╎    ╎    ╎    ╎    ╎    ╎    35     …ase/boot.jl:479; Array
   34╎    ╎    ╎    ╎    ╎    5596   …se/array.jl:844; collect(itr::Base.Genera…
    1╎    ╎    ╎    ╎    ╎     1      …se/array.jl:867; collect_to_with_first!(…
     ╎    ╎    ╎    ╎    ╎     5561   …se/array.jl:870; collect_to_with_first!(…
    2╎    ╎    ╎    ╎    ╎    ╎ 2      …se/array.jl:0; collect_to!(dest::Matrix…
    3╎    ╎    ╎    ╎    ╎    ╎ 3      …se/array.jl:887; collect_to!(dest::Matr…
     ╎    ╎    ╎    ╎    ╎    ╎ 5519   …se/array.jl:892; collect_to!(dest::Matr…
   13╎    ╎    ╎    ╎    ╎    ╎  13     …se/float.jl:0; iterate
     ╎    ╎    ╎    ╎    ╎    ╎  225    …enerator.jl:44; iterate
     ╎    ╎    ╎    ╎    ╎    ╎   225    …terators.jl:1123; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    207    …terators.jl:1110; _piterate1
     ╎    ╎    ╎    ╎    ╎    ╎     7      …se/range.jl:892; iterate
    7╎    ╎    ╎    ╎    ╎    ╎    ╎ 7      @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎    ╎     200    …se/range.jl:894; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 12     …recision.jl:481; unsafe_getindex
   12╎    ╎    ╎    ╎    ╎    ╎    ╎  12     @Base/int.jl:86; -
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 12     …recision.jl:482; unsafe_getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  12     …romotion.jl:423; *
   12╎    ╎    ╎    ╎    ╎    ╎    ╎   12     …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 56     …recision.jl:483; unsafe_getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  37     …recision.jl:84; add12
   19╎    ╎    ╎    ╎    ╎    ╎    ╎   19     …sentials.jl:647; ifelse
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   18     …perators.jl:378; >
   18╎    ╎    ╎    ╎    ╎    ╎    ╎    18     …se/float.jl:536; <
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  19     …recision.jl:85; add12
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   14     …recision.jl:49; canonicalize2
   14╎    ╎    ╎    ╎    ╎    ╎    ╎    14     …se/float.jl:409; +
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   5      …recision.jl:50; canonicalize2
    5╎    ╎    ╎    ╎    ╎    ╎    ╎    5      …se/float.jl:409; +
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 120    …recision.jl:484; unsafe_getindex
  120╎    ╎    ╎    ╎    ╎    ╎    ╎  120    …se/float.jl:409; +
     ╎    ╎    ╎    ╎    ╎    ╎    18     …terators.jl:1114; _piterate1
     ╎    ╎    ╎    ╎    ╎    ╎     18     …terators.jl:1110; _piterate1
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …se/range.jl:892; iterate
    1╎    ╎    ╎    ╎    ╎    ╎    ╎  1      @Base/int.jl:87; +
    3╎    ╎    ╎    ╎    ╎    ╎    ╎ 3      …se/range.jl:893; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 14     …se/range.jl:894; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  3      …recision.jl:482; unsafe_getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   3      …romotion.jl:423; *
    2╎    ╎    ╎    ╎    ╎    ╎    ╎    2      …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …romotion.jl:393; promote
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …romotion.jl:370; _promote
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …e/number.jl:7; convert
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …se/float.jl:159; Float64
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  10     …recision.jl:483; unsafe_getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   6      …recision.jl:84; add12
    3╎    ╎    ╎    ╎    ╎    ╎    ╎    3      …sentials.jl:647; ifelse
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    3      …perators.jl:378; >
    3╎    ╎    ╎    ╎    ╎    ╎    ╎     3      …se/float.jl:536; <
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   4      …recision.jl:85; add12
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …recision.jl:49; canonicalize2
    1╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …se/float.jl:409; +
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    3      …recision.jl:50; canonicalize2
    3╎    ╎    ╎    ╎    ╎    ╎    ╎     3      …se/float.jl:409; +
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …recision.jl:484; unsafe_getindex
    1╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …se/float.jl:409; +
     ╎    ╎    ╎    ╎    ╎    ╎  5281   …enerator.jl:47; iterate
  140╎    ╎    ╎    ╎    ╎    ╎   5281   …pse_comp.jl:62; #284
   11╎    ╎    ╎    ╎    ╎    ╎    11     …ial/trig.jl:98; cos(x::Float64)
    5╎    ╎    ╎    ╎    ╎    ╎    5      …_physics.jl:50; rotation_period(ϕ::F…
  100╎    ╎    ╎    ╎    ╎    ╎    100    …_physics.jl:56; v_scalar(lat::Float6…
 1303╎    ╎    ╎    ╎    ╎    ╎    5025   …_physics.jl:57; v_scalar(lat::Float6…
   44╎    ╎    ╎    ╎    ╎    ╎     44     …se/float.jl:412; /(x::Float64, y::F…
    4╎    ╎    ╎    ╎    ╎    ╎     4      …se/float.jl:0; cos(x::Float64)
 1961╎    ╎    ╎    ╎    ╎    ╎     2037   …perators.jl:587; *
   76╎    ╎    ╎    ╎    ╎    ╎    ╎ 76     …se/float.jl:411; *(x::Float64, y::…
   12╎    ╎    ╎    ╎    ╎    ╎     12     …ial/trig.jl:0; cos(x::Float64)
   16╎    ╎    ╎    ╎    ╎    ╎     16     …ial/trig.jl:98; cos(x::Float64)
     ╎    ╎    ╎    ╎    ╎    ╎     13     …ial/trig.jl:99; cos(x::Float64)
   13╎    ╎    ╎    ╎    ╎    ╎    ╎ 13     …se/float.jl:610; abs
   12╎    ╎    ╎    ╎    ╎    ╎     12     …ial/trig.jl:101; cos(x::Float64)
     ╎    ╎    ╎    ╎    ╎    ╎     151    …ial/trig.jl:104; cos(x::Float64)
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 9      …ial/trig.jl:146; cos_kernel
    9╎    ╎    ╎    ╎    ╎    ╎    ╎  9      …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 53     …ial/trig.jl:147; cos_kernel
   25╎    ╎    ╎    ╎    ╎    ╎    ╎  25     …se/float.jl:409; +
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  19     …ase/math.jl:186; evalpoly
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   19     …ase/math.jl:187; macro expansion
   19╎    ╎    ╎    ╎    ╎    ╎    ╎    19     …se/float.jl:414; muladd
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  9      …perators.jl:587; *
    9╎    ╎    ╎    ╎    ╎    ╎    ╎   9      …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …ial/trig.jl:148; cos_kernel
    1╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 88     …ial/trig.jl:150; cos_kernel
   82╎    ╎    ╎    ╎    ╎    ╎    ╎  82     …se/float.jl:409; +
    6╎    ╎    ╎    ╎    ╎    ╎    ╎  6      …se/float.jl:410; -
     ╎    ╎    ╎    ╎    ╎    ╎     4      …ial/trig.jl:105; cos(x::Float64)
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 4      …se/float.jl:620; isnan
    4╎    ╎    ╎    ╎    ╎    ╎    ╎  4      …se/float.jl:535; !=
     ╎    ╎    ╎    ╎    ╎    ╎     2      …ial/trig.jl:107; cos(x::Float64)
    2╎    ╎    ╎    ╎    ╎    ╎    ╎ 2      …se/float.jl:635; isinf
     ╎    ╎    ╎    ╎    ╎    ╎     12     …ial/trig.jl:110; cos(x::Float64)
    5╎    ╎    ╎    ╎    ╎    ╎    ╎ 5      …rem_pio2.jl:0; rem_pio2_kernel
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 2      …rem_pio2.jl:224; rem_pio2_kernel
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  2      …ase/math.jl:1540; poshighword
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   2      …ase/math.jl:1541; poshighword
    2╎    ╎    ╎    ╎    ╎    ╎    ╎    2      @Base/int.jl:347; &
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 5      …rem_pio2.jl:238; rem_pio2_kernel
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  2      …rem_pio2.jl:52; cody_waite_2c_pio2
    2╎    ╎    ╎    ╎    ╎    ╎    ╎   2      …se/float.jl:414; muladd
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  3      …rem_pio2.jl:54; cody_waite_2c_pio2
    3╎    ╎    ╎    ╎    ╎    ╎    ╎   3      …se/float.jl:414; muladd
    1╎    ╎    ╎    ╎    ╎    ╎     1      …ial/trig.jl:114; cos(x::Float64)
     ╎    ╎    ╎    ╎    ╎    ╎     96     …ial/trig.jl:119; cos(x::Float64)
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 11     …ial/trig.jl:72; sin_kernel
    4╎    ╎    ╎    ╎    ╎    ╎    ╎  4      …se/float.jl:409; +
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  6      …ase/math.jl:186; evalpoly
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   6      …ase/math.jl:187; macro expansion
    6╎    ╎    ╎    ╎    ╎    ╎    ╎    6      …se/float.jl:414; muladd
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …perators.jl:587; *
    1╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 8      …ial/trig.jl:73; sin_kernel
    8╎    ╎    ╎    ╎    ╎    ╎    ╎  8      …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 77     …ial/trig.jl:74; sin_kernel
   77╎    ╎    ╎    ╎    ╎    ╎    ╎  77     …se/float.jl:410; -
   35╎    ╎    ╎    ╎    ╎    ╎     1318   …_physics.jl:50; rotation_period
   23╎    ╎    ╎    ╎    ╎    ╎    ╎ 23     …ase/math.jl:1246; ^(x::Float64, n:…
    4╎    ╎    ╎    ╎    ╎    ╎    ╎ 4      …ial/trig.jl:29; sin(x::Float64)
    4╎    ╎    ╎    ╎    ╎    ╎    ╎ 4      …_physics.jl:0; rotation_period(ϕ::…
   90╎    ╎    ╎    ╎    ╎    ╎    ╎ 90     …_physics.jl:50; rotation_period(ϕ:…
   11╎    ╎    ╎    ╎    ╎    ╎    ╎ 11     …_physics.jl:51; rotation_period(ϕ:…
   16╎    ╎    ╎    ╎    ╎    ╎    ╎ 272    …_physics.jl:52; rotation_period(ϕ:…
    1╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …se/float.jl:0; sin(x::Float64)
    5╎    ╎    ╎    ╎    ╎    ╎    ╎  5      …ial/trig.jl:0; sin(x::Float64)
   23╎    ╎    ╎    ╎    ╎    ╎    ╎  23     …ial/trig.jl:29; sin(x::Float64)
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  7      …ial/trig.jl:30; sin(x::Float64)
    7╎    ╎    ╎    ╎    ╎    ╎    ╎   7      …se/float.jl:610; abs
    5╎    ╎    ╎    ╎    ╎    ╎    ╎  5      …ial/trig.jl:31; sin(x::Float64)
    3╎    ╎    ╎    ╎    ╎    ╎    ╎  3      …ial/trig.jl:32; sin(x::Float64)
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  124    …ial/trig.jl:35; sin(x::Float64)
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   8      …ial/trig.jl:78; sin_kernel
    8╎    ╎    ╎    ╎    ╎    ╎    ╎    8      …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   19     …ial/trig.jl:79; sin_kernel
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …se/float.jl:409; +
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    12     …ase/math.jl:186; evalpoly
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     12     …ase/math.jl:187; macro expansi…
   12╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 12     …se/float.jl:414; muladd
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    6      …perators.jl:587; *
    6╎    ╎    ╎    ╎    ╎    ╎    ╎     6      …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   6      …ial/trig.jl:80; sin_kernel
    6╎    ╎    ╎    ╎    ╎    ╎    ╎    6      …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   91     …ial/trig.jl:81; sin_kernel
   33╎    ╎    ╎    ╎    ╎    ╎    ╎    33     …se/float.jl:411; *
   58╎    ╎    ╎    ╎    ╎    ╎    ╎    58     …se/float.jl:409; +
    2╎    ╎    ╎    ╎    ╎    ╎    ╎  2      …ial/trig.jl:36; sin(x::Float64)
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …ial/trig.jl:38; sin(x::Float64)
    1╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …se/float.jl:635; isinf
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  14     …ial/trig.jl:41; sin(x::Float64)
    3╎    ╎    ╎    ╎    ╎    ╎    ╎   3      …rem_pio2.jl:0; rem_pio2_kernel
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …rem_pio2.jl:224; rem_pio2_kernel
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …ase/math.jl:1540; poshighword
    1╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …sentials.jl:581; reinterpret
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   2      …rem_pio2.jl:226; rem_pio2_kernel
    2╎    ╎    ╎    ╎    ╎    ╎    ╎    2      @Base/int.jl:515; <=
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …rem_pio2.jl:229; rem_pio2_kernel
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   7      …rem_pio2.jl:238; rem_pio2_kernel
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    2      …rem_pio2.jl:52; cody_waite_2c_p…
    2╎    ╎    ╎    ╎    ╎    ╎    ╎     2      …se/float.jl:414; muladd
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    5      …rem_pio2.jl:54; cody_waite_2c_p…
    5╎    ╎    ╎    ╎    ╎    ╎    ╎     5      …se/float.jl:414; muladd
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  2      …ial/trig.jl:48; sin(x::Float64)
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   2      …ial/trig.jl:71; sin_kernel
    2╎    ╎    ╎    ╎    ╎    ╎    ╎    2      …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  69     …ial/trig.jl:50; sin(x::Float64)
    5╎    ╎    ╎    ╎    ╎    ╎    ╎   5      …se/float.jl:407; -
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   10     …ial/trig.jl:139; cos_kernel
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    7      …ase/math.jl:186; evalpoly
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     7      …ase/math.jl:187; macro expansi…
    7╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 7      …se/float.jl:414; muladd
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    2      …perators.jl:587; *
    2╎    ╎    ╎    ╎    ╎    ╎    ╎     2      …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   3      …ial/trig.jl:140; cos_kernel
    3╎    ╎    ╎    ╎    ╎    ╎    ╎    3      …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   2      …ial/trig.jl:141; cos_kernel
    2╎    ╎    ╎    ╎    ╎    ╎    ╎    2      …se/float.jl:410; -
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   49     …ial/trig.jl:142; cos_kernel
   20╎    ╎    ╎    ╎    ╎    ╎    ╎    20     …se/float.jl:411; *
   26╎    ╎    ╎    ╎    ╎    ╎    ╎    26     …se/float.jl:409; +
    3╎    ╎    ╎    ╎    ╎    ╎    ╎    3      …se/float.jl:410; -
   13╎    ╎    ╎    ╎    ╎    ╎    ╎ 879    …_physics.jl:53; rotation_period(ϕ:…
  186╎    ╎    ╎    ╎    ╎    ╎    ╎  186    …se/float.jl:412; /
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  18     …ase/math.jl:1195; ^
   18╎    ╎    ╎    ╎    ╎    ╎    ╎   18     …sentials.jl:581; reinterpret
   42╎    ╎    ╎    ╎    ╎    ╎    ╎  602    …ase/math.jl:1204; ^
   13╎    ╎    ╎    ╎    ╎    ╎    ╎   13     …ase/math.jl:0; ^(x::Float64, n::…
   27╎    ╎    ╎    ╎    ╎    ╎    ╎   27     …ase/math.jl:1246; ^(x::Float64, …
   22╎    ╎    ╎    ╎    ╎    ╎    ╎   479    …ase/math.jl:1248; ^(x::Float64, …
   96╎    ╎    ╎    ╎    ╎    ╎    ╎    96     …ase/math.jl:0; pow_body(x::Floa…
   20╎    ╎    ╎    ╎    ╎    ╎    ╎    20     …ase/math.jl:1251; pow_body(x::F…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    27     …ase/math.jl:1254; pow_body(x::F…
   27╎    ╎    ╎    ╎    ╎    ╎    ╎     27     …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    39     …ase/math.jl:1262; pow_body(x::F…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     39     …perators.jl:378; >
   39╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 39     @Base/int.jl:83; <
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    40     …ase/math.jl:1268; pow_body(x::F…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     40     …perators.jl:587; *
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 40     …romotion.jl:423; *
   40╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  40     …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    36     …ase/math.jl:1270; pow_body(x::F…
   36╎    ╎    ╎    ╎    ╎    ╎    ╎     36     …se/float.jl:409; +
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    26     …ase/math.jl:1271; pow_body(x::F…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     26     @Base/int.jl:538; >>>
   26╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 26     @Base/int.jl:530; >>>
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    14     …ase/math.jl:1273; pow_body(x::F…
   14╎    ╎    ╎    ╎    ╎    ╎    ╎     14     …se/float.jl:414; muladd
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    159    …ase/math.jl:1274; pow_body(x::F…
   86╎    ╎    ╎    ╎    ╎    ╎    ╎     86     …sentials.jl:647; ifelse
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     73     …se/float.jl:623; isfinite
   73╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 73     …se/float.jl:410; -
    1╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …ase/math.jl:0; pow_body(x::Float…
   40╎    ╎    ╎    ╎    ╎    ╎    ╎   40     …ase/math.jl:1251; pow_body(x::Fl…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  60     …perators.jl:587; +
   60╎    ╎    ╎    ╎    ╎    ╎    ╎   60     …se/float.jl:409; +
    6╎    ╎    ╎    ╎    ╎    ╎ 6      …se/array.jl:895; collect_to!(dest::Matr…
     ╎    ╎    ╎    ╎    ╎    ╎ 17     …se/array.jl:896; collect_to!(dest::Matr…
   17╎    ╎    ╎    ╎    ╎    ╎  17     …se/array.jl:1021; setindex!
     ╎    ╎    ╎    ╎    ╎    ╎ 4      …se/array.jl:897; collect_to!(dest::Matr…
    4╎    ╎    ╎    ╎    ╎    ╎  4      @Base/int.jl:87; +
   10╎    ╎    ╎    ╎    ╎    ╎ 10     …_physics.jl:57; v_scalar(lat::Float64, …
    5╎    ╎    ╎    ╎    ╎   5      …se/array.jl:371; _copyto_impl!(dest::Matri…
     ╎    ╎    ╎    ╎    ╎   29     …roadcast.jl:911; materialize!(dest::Matrix…
     ╎    ╎    ╎    ╎    ╎    29     …roadcast.jl:914; materialize!
     ╎    ╎    ╎    ╎    ╎     29     …roadcast.jl:956; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎ 29     …roadcast.jl:997; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎  29     …se/array.jl:388; copyto!
    3╎    ╎    ╎    ╎    ╎    ╎   29     …se/array.jl:368; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎    26     …se/array.jl:376; _copyto_impl!(dest:…
     ╎    ╎    ╎    ╎    ╎    ╎     26     …se/array.jl:337; unsafe_copyto!
   26╎    ╎    ╎    ╎    ╎    ╎    ╎ 26     …ase/cmem.jl:26; memmove
     ╎    ╎    ╎    ╎    ╎  54     …pse_comp.jl:64; eclipse_compute_quantities!…
     ╎    ╎    ╎    ╎    ╎   54     …roadcast.jl:911; materialize!
     ╎    ╎    ╎    ╎    ╎    54     …roadcast.jl:914; materialize!
     ╎    ╎    ╎    ╎    ╎     54     …roadcast.jl:956; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎ 54     …roadcast.jl:1003; copyto!
    2╎    ╎    ╎    ╎    ╎    ╎  3      …simdloop.jl:75; macro expansion
    1╎    ╎    ╎    ╎    ╎    ╎   1      @Base/int.jl:83; <
     ╎    ╎    ╎    ╎    ╎    ╎  49     …simdloop.jl:77; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎   49     …roadcast.jl:1004; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎    13     …roadcast.jl:636; getindex
     ╎    ╎    ╎    ╎    ╎    ╎     9      …roadcast.jl:681; _broadcast_getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 9      …roadcast.jl:705; _getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  9      …roadcast.jl:675; _broadcast_getin…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   9      …ensional.jl:696; getindex
    9╎    ╎    ╎    ╎    ╎    ╎    ╎    9      …sentials.jl:14; getindex
     ╎    ╎    ╎    ╎    ╎    ╎     4      …roadcast.jl:682; _broadcast_getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 4      …roadcast.jl:709; _broadcast_getind…
    4╎    ╎    ╎    ╎    ╎    ╎    ╎  4      …se/float.jl:412; /
     ╎    ╎    ╎    ╎    ╎    ╎    36     …ensional.jl:698; setindex!
   36╎    ╎    ╎    ╎    ╎    ╎     36     …se/array.jl:1024; setindex!
     ╎    ╎    ╎    ╎    ╎    ╎  1      …simdloop.jl:78; macro expansion
    1╎    ╎    ╎    ╎    ╎    ╎   1      @Base/int.jl:87; +
    1╎    ╎    ╎    ╎    ╎    ╎  1      …simdloop.jl:84; macro expansion
   16╎    ╎    ╎    ╎    ╎  4647   …pse_comp.jl:67; eclipse_compute_quantities!…
   10╎    ╎    ╎    ╎    ╎   10     …rraymath.jl:6; -(A::Vector{Float64}, B::Ve…
     ╎    ╎    ╎    ╎    ╎   14     …rraymath.jl:8; -(A::Vector{Float64}, B::Ve…
   14╎    ╎    ╎    ╎    ╎    14     …roadcast.jl:893; broadcast_preserving_zer…
    1╎    ╎    ╎    ╎    ╎   1      …geometry.jl:0; pole_vector_grid!(A::Matrix…
    1╎    ╎    ╎    ╎    ╎   1      …geometry.jl:50; pole_vector_grid!(A::Matri…
   21╎    ╎    ╎    ╎    ╎   4584   …geometry.jl:58; pole_vector_grid!(A::Matri…
   30╎    ╎    ╎    ╎    ╎    30     …se/array.jl:1021; setindex!
     ╎    ╎    ╎    ╎    ╎    1970   …se/array.jl:163; vect
 1970╎    ╎    ╎    ╎    ╎     1970   …ase/boot.jl:477; Array
     ╎    ╎    ╎    ╎    ╎    9      …se/array.jl:165; vect
    9╎    ╎    ╎    ╎    ╎     9      …se/array.jl:1026; __inbounds_setindex!
   68╎    ╎    ╎    ╎    ╎    68     …rraymath.jl:6; -(A::Vector{Float64}, B::V…
     ╎    ╎    ╎    ╎    ╎    12     …rraymath.jl:7; -(A::Vector{Float64}, B::V…
     ╎    ╎    ╎    ╎    ╎     12     …/indices.jl:169; promote_shape
     ╎    ╎    ╎    ╎    ╎    ╎ 12     …/indices.jl:177; promote_shape
     ╎    ╎    ╎    ╎    ╎    ╎  12     …perators.jl:276; !=
     ╎    ╎    ╎    ╎    ╎    ╎   12     …se/range.jl:1134; ==
   12╎    ╎    ╎    ╎    ╎    ╎    12     …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    2364   …rraymath.jl:8; -(A::Vector{Float64}, B::V…
     ╎    ╎    ╎    ╎    ╎     2312   …roadcast.jl:892; broadcast_preserving_ze…
     ╎    ╎    ╎    ╎    ╎    ╎ 2312   …roadcast.jl:903; materialize
     ╎    ╎    ╎    ╎    ╎    ╎  2312   …roadcast.jl:928; copy
     ╎    ╎    ╎    ╎    ╎    ╎   202    …roadcast.jl:956; copyto!
    1╎    ╎    ╎    ╎    ╎    ╎    1      …roadcast.jl:0; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎    15     …roadcast.jl:992; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎     12     …actarray.jl:98; axes
   12╎    ╎    ╎    ╎    ╎    ╎    ╎ 12     …se/array.jl:191; size
     ╎    ╎    ╎    ╎    ╎    ╎     3      …se/tuple.jl:482; ==
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 3      …se/tuple.jl:486; _eq
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  3      …se/range.jl:1134; ==
    3╎    ╎    ╎    ╎    ╎    ╎    ╎   3      …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎    66     …roadcast.jl:1000; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎     66     …roadcast.jl:983; preprocess
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 66     …roadcast.jl:986; preprocess_args
   10╎    ╎    ╎    ╎    ╎    ╎    ╎  10     …se/array.jl:0; size
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  41     …roadcast.jl:984; preprocess
    8╎    ╎    ╎    ╎    ╎    ╎    ╎   25     …roadcast.jl:977; broadcast_unali…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    17     …actarray.jl:1481; unalias
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     17     …actarray.jl:1516; mightalias
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 16     …actarray.jl:1523; _isdisjoint
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  16     …perators.jl:276; !=
   16╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   16     …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …actarray.jl:1540; dataids
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …actarray.jl:1237; pointer
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …/pointer.jl:65; unsafe_conv…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   16     …roadcast.jl:676; extrude
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    16     …roadcast.jl:625; newindexer
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …actarray.jl:98; axes
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …se/array.jl:191; size
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     15     …roadcast.jl:626; shapeindexer
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 15     …roadcast.jl:631; _newindexer
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  15     …perators.jl:276; !=
   15╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   15     …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  15     …roadcast.jl:987; preprocess_args
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   15     …roadcast.jl:984; preprocess
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    15     …roadcast.jl:977; broadcast_unal…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     15     …actarray.jl:1481; unalias
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 15     …actarray.jl:1516; mightalias
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  15     …actarray.jl:1540; dataids
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   15     …actarray.jl:1237; pointer
   15╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    15     …/pointer.jl:65; unsafe_con…
     ╎    ╎    ╎    ╎    ╎    ╎    120    …roadcast.jl:1003; copyto!
    1╎    ╎    ╎    ╎    ╎    ╎     1      …simdloop.jl:0; macro expansion
   65╎    ╎    ╎    ╎    ╎    ╎     65     …simdloop.jl:75; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎     43     …simdloop.jl:77; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 43     …roadcast.jl:1004; macro expansion
   30╎    ╎    ╎    ╎    ╎    ╎    ╎  30     …se/array.jl:1021; setindex!
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  13     …roadcast.jl:636; getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   13     …roadcast.jl:681; _broadcast_geti…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    13     …roadcast.jl:705; _getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     13     …roadcast.jl:675; _broadcast_ge…
   13╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 13     …sentials.jl:13; getindex
     ╎    ╎    ╎    ╎    ╎    ╎     11     …simdloop.jl:78; macro expansion
   11╎    ╎    ╎    ╎    ╎    ╎    ╎ 11     @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎    ╎   2110   …roadcast.jl:223; similar
     ╎    ╎    ╎    ╎    ╎    ╎    2110   …roadcast.jl:224; similar
     ╎    ╎    ╎    ╎    ╎    ╎     2110   …actarray.jl:876; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 2110   …actarray.jl:877; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  2110   …ase/boot.jl:494; Array
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   2110   …ase/boot.jl:486; Array
 2110╎    ╎    ╎    ╎    ╎    ╎    ╎    2110   …ase/boot.jl:477; Array
   38╎    ╎    ╎    ╎    ╎     52     …roadcast.jl:893; broadcast_preserving_ze…
     ╎    ╎    ╎    ╎    ╎    ╎ 14     …roadcast.jl:234; axes
     ╎    ╎    ╎    ╎    ╎    ╎  14     …roadcast.jl:236; _axes
     ╎    ╎    ╎    ╎    ╎    ╎   14     …roadcast.jl:524; combine_axes
     ╎    ╎    ╎    ╎    ╎    ╎    14     …roadcast.jl:543; broadcast_shape
     ╎    ╎    ╎    ╎    ╎    ╎     14     …roadcast.jl:549; _bcs
   14╎    ╎    ╎    ╎    ╎    ╎    ╎ 14     …roadcast.jl:555; _bcs1
  110╎    ╎    ╎    ╎    ╎    110    …sentials.jl:13; getindex
   20╎    ╎    ╎    ╎    ╎   20     …geometry.jl:59; pole_vector_grid!(A::Matri…
    1╎    ╎    ╎    ╎    ╎   1      …geometry.jl:60; pole_vector_grid!(A::Matri…
  129╎    ╎    ╎    ╎    ╎  13278  …pse_comp.jl:70; eclipse_compute_quantities!…
   23╎    ╎    ╎    ╎    ╎   23     …se/array.jl:161; vect(::Float64, ::Vararg{…
    1╎    ╎    ╎    ╎    ╎   1      …se/array.jl:167; vect(::Float64, ::Vararg{…
   11╎    ╎    ╎    ╎    ╎   11     …sentials.jl:0; v_vector(A::Matrix{Vector{F…
   12╎    ╎    ╎    ╎    ╎   12     …geometry.jl:0; v_vector(A::Matrix{Vector{F…
    1╎    ╎    ╎    ╎    ╎   1      …geometry.jl:63; v_vector(A::Matrix{Vector{…
     ╎    ╎    ╎    ╎    ╎   1      …geometry.jl:72; v_vector(A::Matrix{Vector{…
     ╎    ╎    ╎    ╎    ╎    1      …actarray.jl:378; eachindex
     ╎    ╎    ╎    ╎    ╎     1      …actarray.jl:388; eachindex
    1╎    ╎    ╎    ╎    ╎    ╎ 1      …sentials.jl:10; length
  741╎    ╎    ╎    ╎    ╎   5341   …geometry.jl:73; v_vector(A::Matrix{Vector{…
    1╎    ╎    ╎    ╎    ╎    1      …se/array.jl:0; vect(::Float64, ::Vararg{F…
   58╎    ╎    ╎    ╎    ╎    58     …se/array.jl:161; vect(::Float64, ::Vararg…
     ╎    ╎    ╎    ╎    ╎    2055   …se/array.jl:163; vect(::Float64, ::Vararg…
 2043╎    ╎    ╎    ╎    ╎     2043   …ase/boot.jl:477; Array
   12╎    ╎    ╎    ╎    ╎     12     …se/tuple.jl:26; length
   17╎    ╎    ╎    ╎    ╎    17     …se/array.jl:164; vect(::Float64, ::Vararg…
     ╎    ╎    ╎    ╎    ╎    28     …se/array.jl:165; vect(::Float64, ::Vararg…
   16╎    ╎    ╎    ╎    ╎     16     …se/array.jl:1026; __inbounds_setindex!
   12╎    ╎    ╎    ╎    ╎     12     …se/tuple.jl:33; __inbounds_getindex
   14╎    ╎    ╎    ╎    ╎    33     …se/array.jl:166; vect(::Float64, ::Vararg…
     ╎    ╎    ╎    ╎    ╎     19     …se/range.jl:901; iterate
   19╎    ╎    ╎    ╎    ╎    ╎ 19     …romotion.jl:521; ==
   11╎    ╎    ╎    ╎    ╎    11     …se/array.jl:167; vect(::Float64, ::Vararg…
    5╎    ╎    ╎    ╎    ╎    5      …sentials.jl:13; getindex
   67╎    ╎    ╎    ╎    ╎    67     …/generic.jl:310; cross(a::Vector{Float64}…
   15╎    ╎    ╎    ╎    ╎    15     …/generic.jl:311; cross(a::Vector{Float64}…
     ╎    ╎    ╎    ╎    ╎    3      …/generic.jl:314; cross(a::Vector{Float64}…
     ╎    ╎    ╎    ╎    ╎     3      …se/tuple.jl:93; indexed_iterate
     ╎    ╎    ╎    ╎    ╎    ╎ 3      …se/tuple.jl:93; indexed_iterate
    3╎    ╎    ╎    ╎    ╎    ╎  3      …sentials.jl:13; getindex
     ╎    ╎    ╎    ╎    ╎    33     …/generic.jl:315; cross(a::Vector{Float64}…
     ╎    ╎    ╎    ╎    ╎     33     …se/tuple.jl:93; indexed_iterate
   14╎    ╎    ╎    ╎    ╎    ╎ 14     …sentials.jl:13; getindex
     ╎    ╎    ╎    ╎    ╎    ╎ 19     …se/tuple.jl:93; indexed_iterate
   19╎    ╎    ╎    ╎    ╎    ╎  19     …sentials.jl:13; getindex
     ╎    ╎    ╎    ╎    ╎    2274   …/generic.jl:316; cross(a::Vector{Float64}…
     ╎    ╎    ╎    ╎    ╎     2164   …se/array.jl:163; vect
 2164╎    ╎    ╎    ╎    ╎    ╎ 2164   …ase/boot.jl:477; Array
     ╎    ╎    ╎    ╎    ╎     11     …se/array.jl:165; vect
   11╎    ╎    ╎    ╎    ╎    ╎ 11     …se/array.jl:1026; __inbounds_setindex!
   22╎    ╎    ╎    ╎    ╎     22     …se/array.jl:167; vect
   48╎    ╎    ╎    ╎    ╎     48     …se/float.jl:411; *
   29╎    ╎    ╎    ╎    ╎     29     …se/float.jl:410; -
 2437╎    ╎    ╎    ╎    ╎   3414   …geometry.jl:74; v_vector(A::Matrix{Vector{…
   21╎    ╎    ╎    ╎    ╎    21     …roadcast.jl:1343; broadcasted(::typeof(/)…
     ╎    ╎    ╎    ╎    ╎    25     …roadcast.jl:1347; broadcasted(::typeof(/)…
     ╎    ╎    ╎    ╎    ╎     25     …roadcast.jl:1349; broadcasted
     ╎    ╎    ╎    ╎    ╎    ╎ 25     …roadcast.jl:178; Broadcasted
   25╎    ╎    ╎    ╎    ╎    ╎  25     …roadcast.jl:178; Broadcasted
   69╎    ╎    ╎    ╎    ╎    69     …roadcast.jl:910; materialize!(dest::Vecto…
   24╎    ╎    ╎    ╎    ╎    229    …roadcast.jl:911; materialize!(dest::Vecto…
     ╎    ╎    ╎    ╎    ╎     205    …roadcast.jl:914; materialize!
     ╎    ╎    ╎    ╎    ╎    ╎ 6      …actarray.jl:98; axes
    6╎    ╎    ╎    ╎    ╎    ╎  6      …se/array.jl:191; size
     ╎    ╎    ╎    ╎    ╎    ╎ 189    …roadcast.jl:956; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎  22     …roadcast.jl:1000; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎   22     …roadcast.jl:983; preprocess
     ╎    ╎    ╎    ╎    ╎    ╎    22     …roadcast.jl:986; preprocess_args
     ╎    ╎    ╎    ╎    ╎    ╎     22     …roadcast.jl:984; preprocess
    6╎    ╎    ╎    ╎    ╎    ╎    ╎ 6      …roadcast.jl:977; broadcast_unalias
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 16     …roadcast.jl:676; extrude
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  16     …roadcast.jl:625; newindexer
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   16     …roadcast.jl:626; shapeindexer
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    16     …roadcast.jl:631; _newindexer
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     16     …perators.jl:276; !=
   16╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 16     …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎  167    …roadcast.jl:1003; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎   13     …simdloop.jl:72; macro expansion
   13╎    ╎    ╎    ╎    ╎    ╎    13     @Base/int.jl:83; <
   20╎    ╎    ╎    ╎    ╎    ╎   20     …simdloop.jl:75; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎   110    …simdloop.jl:77; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎    110    …roadcast.jl:1004; macro expansion
   99╎    ╎    ╎    ╎    ╎    ╎     99     …se/array.jl:1021; setindex!
     ╎    ╎    ╎    ╎    ╎    ╎     11     …roadcast.jl:636; getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 6      …roadcast.jl:681; _broadcast_getind…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  6      …roadcast.jl:705; _getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   6      …roadcast.jl:675; _broadcast_geti…
    6╎    ╎    ╎    ╎    ╎    ╎    ╎    6      …sentials.jl:13; getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 5      …roadcast.jl:682; _broadcast_getind…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  5      …roadcast.jl:709; _broadcast_getin…
    5╎    ╎    ╎    ╎    ╎    ╎    ╎   5      …se/float.jl:412; /
     ╎    ╎    ╎    ╎    ╎    ╎   24     …simdloop.jl:78; macro expansion
   24╎    ╎    ╎    ╎    ╎    ╎    24     @Base/int.jl:87; +
    8╎    ╎    ╎    ╎    ╎    ╎ 10     …roadcast.jl:309; instantiate
     ╎    ╎    ╎    ╎    ╎    ╎  2      …roadcast.jl:585; check_broadcast_axes
     ╎    ╎    ╎    ╎    ╎    ╎   2      …roadcast.jl:582; check_broadcast_axes
     ╎    ╎    ╎    ╎    ╎    ╎    2      …actarray.jl:98; axes
    2╎    ╎    ╎    ╎    ╎    ╎     2      …se/array.jl:191; size
    8╎    ╎    ╎    ╎    ╎    8      @Base/int.jl:0; materialize!(dest::Vector{…
   25╎    ╎    ╎    ╎    ╎    25     …/generic.jl:595; norm(itr::Vector{Float64…
   31╎    ╎    ╎    ╎    ╎    600    …/generic.jl:596; norm(itr::Vector{Float64…
   11╎    ╎    ╎    ╎    ╎     11     …/generic.jl:0; generic_norm2(x::Vector{F…
    1╎    ╎    ╎    ╎    ╎     1      …/generic.jl:462; generic_norm2(x::Vector…
   12╎    ╎    ╎    ╎    ╎     12     …/generic.jl:595; norm(itr::Vector{Float6…
     ╎    ╎    ╎    ╎    ╎     7      …/generic.jl:596; norm(itr::Vector{Float6…
     ╎    ╎    ╎    ╎    ╎    ╎ 7      …actarray.jl:1220; isempty
    7╎    ╎    ╎    ╎    ╎    ╎  7      …romotion.jl:521; ==
    3╎    ╎    ╎    ╎    ╎     3      …/generic.jl:597; norm(itr::Vector{Float6…
     ╎    ╎    ╎    ╎    ╎     535    …/generic.jl:598; norm(itr::Vector{Float6…
    2╎    ╎    ╎    ╎    ╎    ╎ 535    …rc/dense.jl:106; norm2
   21╎    ╎    ╎    ╎    ╎    ╎  21     …se/array.jl:0; generic_norm2(x::Vector…
    2╎    ╎    ╎    ╎    ╎    ╎  2      @Base/int.jl:83; <
    1╎    ╎    ╎    ╎    ╎    ╎  1      …e/reduce.jl:0; _mapreduce(f::typeof(Li…
   14╎    ╎    ╎    ╎    ╎    ╎  14     …e/reduce.jl:428; _mapreduce(f::typeof(…
   33╎    ╎    ╎    ╎    ╎    ╎  33     …/generic.jl:0; generic_norm2(x::Vector…
   16╎    ╎    ╎    ╎    ╎    ╎  16     …/generic.jl:462; generic_norm2(x::Vect…
     ╎    ╎    ╎    ╎    ╎    ╎  228    …/generic.jl:463; generic_norm2(x::Vect…
     ╎    ╎    ╎    ╎    ╎    ╎   228    …/generic.jl:527; normInf
     ╎    ╎    ╎    ╎    ╎    ╎    228    …/generic.jl:453; generic_normInf
     ╎    ╎    ╎    ╎    ╎    ╎     228    …educedim.jl:357; mapreduce
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 228    …educedim.jl:357; #mapreduce#821
   15╎    ╎    ╎    ╎    ╎    ╎    ╎  228    …educedim.jl:365; _mapreduce_dim
   21╎    ╎    ╎    ╎    ╎    ╎    ╎   21     …e/reduce.jl:0; _mapreduce(f::typ…
   16╎    ╎    ╎    ╎    ╎    ╎    ╎   16     …e/reduce.jl:428; _mapreduce(f::t…
   10╎    ╎    ╎    ╎    ╎    ╎    ╎   10     …e/reduce.jl:431; _mapreduce(f::t…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   10     …e/reduce.jl:438; _mapreduce(f::t…
   10╎    ╎    ╎    ╎    ╎    ╎    ╎    10     …sentials.jl:13; getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   41     …e/reduce.jl:440; _mapreduce(f::t…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    3      …ase/math.jl:908; max
    3╎    ╎    ╎    ╎    ╎    ╎    ╎     3      …se/float.jl:410; -
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    34     …ase/math.jl:909; max
   15╎    ╎    ╎    ╎    ╎    ╎    ╎     15     …sentials.jl:647; ifelse
   19╎    ╎    ╎    ╎    ╎    ╎    ╎     19     …oatfuncs.jl:15; signbit
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    4      …ase/math.jl:911; max
    4╎    ╎    ╎    ╎    ╎    ╎    ╎     4      …sentials.jl:647; ifelse
   31╎    ╎    ╎    ╎    ╎    ╎    ╎   41     …e/reduce.jl:441; _mapreduce(f::t…
   10╎    ╎    ╎    ╎    ╎    ╎    ╎    10     @Base/int.jl:83; <
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   74     …e/reduce.jl:443; _mapreduce(f::t…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    14     …ase/math.jl:908; max
   14╎    ╎    ╎    ╎    ╎    ╎    ╎     14     …se/float.jl:410; -
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    44     …ase/math.jl:909; max
   30╎    ╎    ╎    ╎    ╎    ╎    ╎     30     …sentials.jl:647; ifelse
   14╎    ╎    ╎    ╎    ╎    ╎    ╎     14     …oatfuncs.jl:15; signbit
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    16     …ase/math.jl:911; max
   16╎    ╎    ╎    ╎    ╎    ╎    ╎     16     …sentials.jl:647; ifelse
    9╎    ╎    ╎    ╎    ╎    ╎  34     …/generic.jl:464; generic_norm2(x::Vect…
   14╎    ╎    ╎    ╎    ╎    ╎   25     …se/float.jl:635; isinf
   11╎    ╎    ╎    ╎    ╎    ╎    11     …se/float.jl:610; abs
     ╎    ╎    ╎    ╎    ╎    ╎  13     …/generic.jl:465; generic_norm2(x::Vect…
     ╎    ╎    ╎    ╎    ╎    ╎   13     …se/array.jl:945; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    13     …se/array.jl:945; iterate
   13╎    ╎    ╎    ╎    ╎    ╎     13     …sentials.jl:13; getindex
   18╎    ╎    ╎    ╎    ╎    ╎  110    …/generic.jl:467; generic_norm2(x::Vect…
     ╎    ╎    ╎    ╎    ╎    ╎   75     …se/float.jl:623; isfinite
   68╎    ╎    ╎    ╎    ╎    ╎    68     …se/float.jl:410; -
     ╎    ╎    ╎    ╎    ╎    ╎    7      …se/float.jl:620; isnan
    7╎    ╎    ╎    ╎    ╎    ╎     7      …se/float.jl:535; !=
     ╎    ╎    ╎    ╎    ╎    ╎   17     …perators.jl:587; *
   17╎    ╎    ╎    ╎    ╎    ╎    17     …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎  11     …/generic.jl:468; generic_norm2(x::Vect…
     ╎    ╎    ╎    ╎    ╎    ╎   11     …/generic.jl:459; norm_sqr
     ╎    ╎    ╎    ╎    ╎    ╎    11     …e/number.jl:189; abs2
   11╎    ╎    ╎    ╎    ╎    ╎     11     …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎  44     …/generic.jl:470; generic_norm2(x::Vect…
   15╎    ╎    ╎    ╎    ╎    ╎   44     …se/array.jl:945; iterate
   14╎    ╎    ╎    ╎    ╎    ╎    14     …sentials.jl:13; getindex
   15╎    ╎    ╎    ╎    ╎    ╎    15     @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎    ╎  6      …/generic.jl:473; generic_norm2(x::Vect…
    6╎    ╎    ╎    ╎    ╎    ╎   6      …se/float.jl:409; +
 2049╎    ╎    ╎    ╎    ╎   2365   …geometry.jl:75; v_vector(A::Matrix{Vector{…
   24╎    ╎    ╎    ╎    ╎    24     …roadcast.jl:1343; broadcasted(::typeof(*)…
     ╎    ╎    ╎    ╎    ╎    23     …roadcast.jl:1347; broadcasted(::typeof(*)…
     ╎    ╎    ╎    ╎    ╎     23     …roadcast.jl:1349; broadcasted
     ╎    ╎    ╎    ╎    ╎    ╎ 23     …roadcast.jl:178; Broadcasted
   23╎    ╎    ╎    ╎    ╎    ╎  23     …roadcast.jl:178; Broadcasted
   67╎    ╎    ╎    ╎    ╎    67     …roadcast.jl:910; materialize!(dest::Vecto…
   41╎    ╎    ╎    ╎    ╎    179    …roadcast.jl:911; materialize!(dest::Vecto…
     ╎    ╎    ╎    ╎    ╎     138    …roadcast.jl:914; materialize!
     ╎    ╎    ╎    ╎    ╎    ╎ 7      …actarray.jl:98; axes
    7╎    ╎    ╎    ╎    ╎    ╎  7      …se/array.jl:191; size
     ╎    ╎    ╎    ╎    ╎    ╎ 124    …roadcast.jl:956; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎  16     …roadcast.jl:1000; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎   16     …roadcast.jl:983; preprocess
     ╎    ╎    ╎    ╎    ╎    ╎    16     …roadcast.jl:986; preprocess_args
     ╎    ╎    ╎    ╎    ╎    ╎     16     …roadcast.jl:984; preprocess
    3╎    ╎    ╎    ╎    ╎    ╎    ╎ 3      …roadcast.jl:977; broadcast_unalias
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 13     …roadcast.jl:676; extrude
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  13     …roadcast.jl:625; newindexer
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   13     …roadcast.jl:626; shapeindexer
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    13     …roadcast.jl:631; _newindexer
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     13     …perators.jl:276; !=
   13╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 13     …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎  108    …roadcast.jl:1003; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎   14     …simdloop.jl:72; macro expansion
   14╎    ╎    ╎    ╎    ╎    ╎    14     @Base/int.jl:83; <
   47╎    ╎    ╎    ╎    ╎    ╎   47     …simdloop.jl:75; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎   45     …simdloop.jl:77; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎    45     …roadcast.jl:1004; macro expansion
   33╎    ╎    ╎    ╎    ╎    ╎     33     …se/array.jl:1021; setindex!
     ╎    ╎    ╎    ╎    ╎    ╎     12     …roadcast.jl:636; getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 12     …roadcast.jl:682; _broadcast_getind…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  12     …roadcast.jl:709; _broadcast_getin…
   12╎    ╎    ╎    ╎    ╎    ╎    ╎   12     …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎   2      …simdloop.jl:78; macro expansion
    2╎    ╎    ╎    ╎    ╎    ╎    2      @Base/int.jl:87; +
    7╎    ╎    ╎    ╎    ╎    ╎ 7      …roadcast.jl:309; instantiate
   20╎    ╎    ╎    ╎    ╎    20     …sentials.jl:13; getindex
    3╎    ╎    ╎    ╎    ╎    3      @Base/int.jl:0; materialize!(dest::Vector{…
 1893╎    ╎    ╎    ╎    ╎   1979   …geometry.jl:76; v_vector(A::Matrix{Vector{…
   86╎    ╎    ╎    ╎    ╎    86     …se/array.jl:1021; setindex!(A::Matrix{Vec…
    1╎    ╎    ╎    ╎    ╎   1      …geometry.jl:78; v_vector(A::Matrix{Vector{…
     ╎    ╎    ╎    ╎    ╎  4343   …pse_comp.jl:73; eclipse_compute_quantities!…
     ╎    ╎    ╎    ╎    ╎   17     …roadcast.jl:1244; dotview
     ╎    ╎    ╎    ╎    ╎    17     …se/views.jl:149; maybeview
   17╎    ╎    ╎    ╎    ╎     17     …sentials.jl:13; getindex
     ╎    ╎    ╎    ╎    ╎   332    …roadcast.jl:911; materialize!
     ╎    ╎    ╎    ╎    ╎    332    …roadcast.jl:914; materialize!
     ╎    ╎    ╎    ╎    ╎     3      …actarray.jl:98; axes
    3╎    ╎    ╎    ╎    ╎    ╎ 3      …se/array.jl:191; size
     ╎    ╎    ╎    ╎    ╎     325    …roadcast.jl:956; copyto!
    3╎    ╎    ╎    ╎    ╎    ╎ 3      …roadcast.jl:0; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎ 7      …roadcast.jl:996; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎  7      …se/tuple.jl:482; ==
     ╎    ╎    ╎    ╎    ╎    ╎   7      …se/tuple.jl:486; _eq
     ╎    ╎    ╎    ╎    ╎    ╎    7      …se/range.jl:1134; ==
    7╎    ╎    ╎    ╎    ╎    ╎     7      …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎ 315    …roadcast.jl:997; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎  315    …se/array.jl:388; copyto!
   45╎    ╎    ╎    ╎    ╎    ╎   315    …se/array.jl:368; copyto!
   12╎    ╎    ╎    ╎    ╎    ╎    12     …se/array.jl:0; _copyto_impl!(dest::V…
   44╎    ╎    ╎    ╎    ╎    ╎    44     …se/array.jl:371; _copyto_impl!(dest:…
     ╎    ╎    ╎    ╎    ╎    ╎    46     …se/array.jl:372; _copyto_impl!(dest:…
   46╎    ╎    ╎    ╎    ╎    ╎     46     …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎    47     …se/array.jl:374; _copyto_impl!(dest:…
    9╎    ╎    ╎    ╎    ╎    ╎     9      …actarray.jl:0; checkbounds
   12╎    ╎    ╎    ╎    ╎    ╎     12     …actarray.jl:700; checkbounds
     ╎    ╎    ╎    ╎    ╎    ╎     13     …actarray.jl:702; checkbounds
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 13     …actarray.jl:687; checkbounds
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  13     …actarray.jl:768; checkindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   13     …actarray.jl:763; checkindex
   12╎    ╎    ╎    ╎    ╎    ╎    ╎    12     @Base/int.jl:86; -
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    1      @Base/int.jl:513; <
    7╎    ╎    ╎    ╎    ╎    ╎     7      @Base/int.jl:86; -
     ╎    ╎    ╎    ╎    ╎    ╎     6      …se/range.jl:5; Colon
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 6      …se/range.jl:403; UnitRange
    6╎    ╎    ╎    ╎    ╎    ╎    ╎  6      …se/range.jl:414; unitrange_last
     ╎    ╎    ╎    ╎    ╎    ╎    37     …se/array.jl:375; _copyto_impl!(dest:…
    7╎    ╎    ╎    ╎    ╎    ╎     7      …actarray.jl:700; checkbounds
     ╎    ╎    ╎    ╎    ╎    ╎     24     …actarray.jl:702; checkbounds
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 24     …actarray.jl:687; checkbounds
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  11     …actarray.jl:768; checkindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   11     …actarray.jl:763; checkindex
   11╎    ╎    ╎    ╎    ╎    ╎    ╎    11     @Base/int.jl:513; <
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  13     …actarray.jl:389; eachindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   13     …actarray.jl:137; axes1
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    13     …actarray.jl:98; axes
   13╎    ╎    ╎    ╎    ╎    ╎    ╎     13     …se/array.jl:191; size
     ╎    ╎    ╎    ╎    ╎    ╎     6      …se/range.jl:5; Colon
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 6      …se/range.jl:403; UnitRange
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  6      …se/range.jl:414; unitrange_last
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   6      …perators.jl:425; >=
    6╎    ╎    ╎    ╎    ╎    ╎    ╎    6      @Base/int.jl:514; <=
     ╎    ╎    ╎    ╎    ╎    ╎    84     …se/array.jl:376; _copyto_impl!(dest:…
     ╎    ╎    ╎    ╎    ╎    ╎     9      …se/array.jl:331; unsafe_copyto!
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 9      …actarray.jl:1240; pointer
    9╎    ╎    ╎    ╎    ╎    ╎    ╎  9      …/pointer.jl:282; +
     ╎    ╎    ╎    ╎    ╎    ╎     14     …se/array.jl:332; unsafe_copyto!
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 14     …actarray.jl:1240; pointer
   14╎    ╎    ╎    ╎    ╎    ╎    ╎  14     …/pointer.jl:282; +
     ╎    ╎    ╎    ╎    ╎    ╎     61     …se/array.jl:337; unsafe_copyto!
   47╎    ╎    ╎    ╎    ╎    ╎    ╎ 50     …ase/cmem.jl:26; memmove
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  3      …sentials.jl:543; cconvert
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   3      …e/number.jl:7; convert
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    3      …ase/boot.jl:789; UInt64
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     3      …ase/boot.jl:759; toUInt64
    3╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 3      …ase/boot.jl:648; check_top_bit
   11╎    ╎    ╎    ╎    ╎    ╎    ╎ 11     @Base/int.jl:88; *
     ╎    ╎    ╎    ╎    ╎     4      …roadcast.jl:309; instantiate
     ╎    ╎    ╎    ╎    ╎    ╎ 4      …roadcast.jl:582; check_broadcast_axes
     ╎    ╎    ╎    ╎    ╎    ╎  4      …actarray.jl:98; axes
    4╎    ╎    ╎    ╎    ╎    ╎   4      …se/array.jl:191; size
    8╎    ╎    ╎    ╎    ╎   8      …sentials.jl:13; getindex
     ╎    ╎    ╎    ╎    ╎   3986   …c/matmul.jl:53; *
     ╎    ╎    ╎    ╎    ╎    2134   …actarray.jl:833; similar
     ╎    ╎    ╎    ╎    ╎     2134   …se/array.jl:420; similar
     ╎    ╎    ╎    ╎    ╎    ╎ 2134   …ase/boot.jl:486; Array
 2134╎    ╎    ╎    ╎    ╎    ╎  2134   …ase/boot.jl:477; Array
   16╎    ╎    ╎    ╎    ╎    16     …se/array.jl:190; size
     ╎    ╎    ╎    ╎    ╎    1836   …c/matmul.jl:237; mul!
     ╎    ╎    ╎    ╎    ╎     1836   …c/matmul.jl:66; mul!
   26╎    ╎    ╎    ╎    ╎    ╎ 1836   …c/matmul.jl:71; generic_matvecmul!
    1╎    ╎    ╎    ╎    ╎    ╎  1      …src/blas.jl:643; gemv!(trans::Char, al…
    1╎    ╎    ╎    ╎    ╎    ╎  1      …src/blas.jl:667; gemv!(trans::Char, al…
   70╎    ╎    ╎    ╎    ╎    ╎  70     …c/matmul.jl:0; gemv!(y::Vector{Float64…
   58╎    ╎    ╎    ╎    ╎    ╎  58     …c/matmul.jl:401; gemv!(y::Vector{Float…
     ╎    ╎    ╎    ╎    ╎    ╎  14     …c/matmul.jl:403; gemv!(y::Vector{Float…
     ╎    ╎    ╎    ╎    ╎    ╎   14     …c/matmul.jl:656; lapack_size
     ╎    ╎    ╎    ╎    ╎    ╎    14     …ase/char.jl:213; ==
   14╎    ╎    ╎    ╎    ╎    ╎     14     …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎  21     …c/matmul.jl:404; gemv!(y::Vector{Float…
    6╎    ╎    ╎    ╎    ╎    ╎   6      …sentials.jl:10; length
     ╎    ╎    ╎    ╎    ╎    ╎   15     …perators.jl:276; !=
   15╎    ╎    ╎    ╎    ╎    ╎    15     …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎  10     …c/matmul.jl:409; gemv!(y::Vector{Float…
   10╎    ╎    ╎    ╎    ╎    ╎   10     …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎  13     …c/matmul.jl:410; gemv!(y::Vector{Float…
     ╎    ╎    ╎    ╎    ╎    ╎   13     …romotion.jl:399; promote
     ╎    ╎    ╎    ╎    ╎    ╎    13     …romotion.jl:377; _promote
     ╎    ╎    ╎    ╎    ╎    ╎     13     …e/number.jl:7; convert
   13╎    ╎    ╎    ╎    ╎    ╎    ╎ 13     …se/float.jl:165; Float64
     ╎    ╎    ╎    ╎    ╎    ╎  10     …c/matmul.jl:411; gemv!(y::Vector{Float…
   10╎    ╎    ╎    ╎    ╎    ╎   10     …perators.jl:1296; in
   55╎    ╎    ╎    ╎    ╎    ╎  1612   …c/matmul.jl:415; gemv!(y::Vector{Float…
   70╎    ╎    ╎    ╎    ╎    ╎   70     …src/blas.jl:643; gemv!(trans::Char, a…
     ╎    ╎    ╎    ╎    ╎    ╎   13     …src/blas.jl:647; gemv!(trans::Char, a…
   13╎    ╎    ╎    ╎    ╎    ╎    13     …se/array.jl:190; size
   25╎    ╎    ╎    ╎    ╎    ╎   25     …src/blas.jl:648; gemv!(trans::Char, a…
     ╎    ╎    ╎    ╎    ╎    ╎   1      …src/blas.jl:657; gemv!(trans::Char, a…
     ╎    ╎    ╎    ╎    ╎    ╎    1      …src/blas.jl:177; vec_pointer_stride
     ╎    ╎    ╎    ╎    ╎    ╎     1      …actarray.jl:1237; pointer
    1╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …/pointer.jl:65; unsafe_convert
     ╎    ╎    ╎    ╎    ╎    ╎   15     …src/blas.jl:659; gemv!(trans::Char, a…
     ╎    ╎    ╎    ╎    ╎    ╎    15     …actarray.jl:1237; pointer
   15╎    ╎    ╎    ╎    ╎    ╎     15     …/pointer.jl:65; unsafe_convert
     ╎    ╎    ╎    ╎    ╎    ╎   31     …src/blas.jl:666; gemv!(trans::Char, a…
     ╎    ╎    ╎    ╎    ╎    ╎    31     …perators.jl:587; max
     ╎    ╎    ╎    ╎    ╎    ╎     31     …romotion.jl:532; max
   18╎    ╎    ╎    ╎    ╎    ╎    ╎ 18     …sentials.jl:647; ifelse
   13╎    ╎    ╎    ╎    ╎    ╎    ╎ 13     @Base/int.jl:83; <
 1380╎    ╎    ╎    ╎    ╎    ╎   1402   …src/blas.jl:667; gemv!(trans::Char, a…
     ╎    ╎    ╎    ╎    ╎    ╎    22     …sentials.jl:543; cconvert
     ╎    ╎    ╎    ╎    ╎    ╎     22     …fpointer.jl:105; convert
   22╎    ╎    ╎    ╎    ╎    ╎    ╎ 22     …refvalue.jl:8; RefValue
     ╎    ╎    ╎    ╎    ╎  4479   …pse_comp.jl:74; eclipse_compute_quantities!…
     ╎    ╎    ╎    ╎    ╎   56     …roadcast.jl:1244; dotview
     ╎    ╎    ╎    ╎    ╎    56     …se/views.jl:149; maybeview
   56╎    ╎    ╎    ╎    ╎     56     …sentials.jl:13; getindex
     ╎    ╎    ╎    ╎    ╎   334    …roadcast.jl:911; materialize!
     ╎    ╎    ╎    ╎    ╎    334    …roadcast.jl:914; materialize!
     ╎    ╎    ╎    ╎    ╎     7      …actarray.jl:98; axes
    7╎    ╎    ╎    ╎    ╎    ╎ 7      …se/array.jl:191; size
     ╎    ╎    ╎    ╎    ╎     319    …roadcast.jl:956; copyto!
    1╎    ╎    ╎    ╎    ╎    ╎ 1      …roadcast.jl:0; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎ 8      …roadcast.jl:996; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎  8      …se/tuple.jl:482; ==
     ╎    ╎    ╎    ╎    ╎    ╎   8      …se/tuple.jl:486; _eq
     ╎    ╎    ╎    ╎    ╎    ╎    8      …se/range.jl:1134; ==
    8╎    ╎    ╎    ╎    ╎    ╎     8      …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎ 310    …roadcast.jl:997; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎  310    …se/array.jl:388; copyto!
   50╎    ╎    ╎    ╎    ╎    ╎   310    …se/array.jl:368; copyto!
    8╎    ╎    ╎    ╎    ╎    ╎    8      …se/array.jl:0; _copyto_impl!(dest::V…
   41╎    ╎    ╎    ╎    ╎    ╎    41     …se/array.jl:371; _copyto_impl!(dest:…
     ╎    ╎    ╎    ╎    ╎    ╎    39     …se/array.jl:372; _copyto_impl!(dest:…
   39╎    ╎    ╎    ╎    ╎    ╎     39     …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎    40     …se/array.jl:374; _copyto_impl!(dest:…
    9╎    ╎    ╎    ╎    ╎    ╎     9      …actarray.jl:0; checkbounds
   10╎    ╎    ╎    ╎    ╎    ╎     10     …actarray.jl:700; checkbounds
     ╎    ╎    ╎    ╎    ╎    ╎     5      …actarray.jl:702; checkbounds
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 5      …actarray.jl:687; checkbounds
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  5      …actarray.jl:768; checkindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   5      …actarray.jl:763; checkindex
    5╎    ╎    ╎    ╎    ╎    ╎    ╎    5      @Base/int.jl:86; -
    1╎    ╎    ╎    ╎    ╎    ╎     1      @Base/int.jl:87; +
   15╎    ╎    ╎    ╎    ╎    ╎     15     @Base/int.jl:86; -
     ╎    ╎    ╎    ╎    ╎    ╎    51     …se/array.jl:375; _copyto_impl!(dest:…
    5╎    ╎    ╎    ╎    ╎    ╎     5      …actarray.jl:700; checkbounds
     ╎    ╎    ╎    ╎    ╎    ╎     34     …actarray.jl:702; checkbounds
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 34     …actarray.jl:687; checkbounds
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  20     …actarray.jl:768; checkindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   20     …actarray.jl:763; checkindex
   20╎    ╎    ╎    ╎    ╎    ╎    ╎    20     @Base/int.jl:513; <
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  14     …actarray.jl:389; eachindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   14     …actarray.jl:137; axes1
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    14     …actarray.jl:98; axes
   14╎    ╎    ╎    ╎    ╎    ╎    ╎     14     …se/array.jl:191; size
     ╎    ╎    ╎    ╎    ╎    ╎     12     …se/range.jl:5; Colon
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 12     …se/range.jl:403; UnitRange
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  12     …se/range.jl:414; unitrange_last
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   12     …perators.jl:425; >=
   12╎    ╎    ╎    ╎    ╎    ╎    ╎    12     @Base/int.jl:514; <=
     ╎    ╎    ╎    ╎    ╎    ╎    81     …se/array.jl:376; _copyto_impl!(dest:…
     ╎    ╎    ╎    ╎    ╎    ╎     12     …se/array.jl:331; unsafe_copyto!
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 12     …actarray.jl:1240; pointer
   10╎    ╎    ╎    ╎    ╎    ╎    ╎  10     …/pointer.jl:282; +
    2╎    ╎    ╎    ╎    ╎    ╎    ╎  2      …/pointer.jl:65; unsafe_convert
     ╎    ╎    ╎    ╎    ╎    ╎     10     …se/array.jl:332; unsafe_copyto!
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 10     …actarray.jl:1240; pointer
   10╎    ╎    ╎    ╎    ╎    ╎    ╎  10     …/pointer.jl:282; +
     ╎    ╎    ╎    ╎    ╎    ╎     59     …se/array.jl:337; unsafe_copyto!
   47╎    ╎    ╎    ╎    ╎    ╎    ╎ 48     …ase/cmem.jl:26; memmove
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …sentials.jl:543; cconvert
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …e/number.jl:7; convert
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …ase/boot.jl:789; UInt64
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …ase/boot.jl:759; toUInt64
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …ase/boot.jl:648; check_top_bit
   11╎    ╎    ╎    ╎    ╎    ╎    ╎ 11     @Base/int.jl:88; *
     ╎    ╎    ╎    ╎    ╎     8      …roadcast.jl:309; instantiate
     ╎    ╎    ╎    ╎    ╎    ╎ 8      …roadcast.jl:582; check_broadcast_axes
     ╎    ╎    ╎    ╎    ╎    ╎  8      …actarray.jl:98; axes
    8╎    ╎    ╎    ╎    ╎    ╎   8      …se/array.jl:191; size
   14╎    ╎    ╎    ╎    ╎   14     …sentials.jl:13; getindex
     ╎    ╎    ╎    ╎    ╎   4075   …c/matmul.jl:53; *
     ╎    ╎    ╎    ╎    ╎    2056   …actarray.jl:833; similar
     ╎    ╎    ╎    ╎    ╎     2056   …se/array.jl:420; similar
     ╎    ╎    ╎    ╎    ╎    ╎ 2056   …ase/boot.jl:486; Array
 2056╎    ╎    ╎    ╎    ╎    ╎  2056   …ase/boot.jl:477; Array
     ╎    ╎    ╎    ╎    ╎    2019   …c/matmul.jl:237; mul!
     ╎    ╎    ╎    ╎    ╎     2019   …c/matmul.jl:66; mul!
   31╎    ╎    ╎    ╎    ╎    ╎ 2019   …c/matmul.jl:71; generic_matvecmul!
    3╎    ╎    ╎    ╎    ╎    ╎  3      …src/blas.jl:643; gemv!(trans::Char, al…
   65╎    ╎    ╎    ╎    ╎    ╎  65     …c/matmul.jl:0; gemv!(y::Vector{Float64…
   57╎    ╎    ╎    ╎    ╎    ╎  57     …c/matmul.jl:401; gemv!(y::Vector{Float…
     ╎    ╎    ╎    ╎    ╎    ╎  19     …c/matmul.jl:403; gemv!(y::Vector{Float…
     ╎    ╎    ╎    ╎    ╎    ╎   19     …c/matmul.jl:656; lapack_size
     ╎    ╎    ╎    ╎    ╎    ╎    19     …ase/char.jl:213; ==
   19╎    ╎    ╎    ╎    ╎    ╎     19     …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎  67     …c/matmul.jl:404; gemv!(y::Vector{Float…
   22╎    ╎    ╎    ╎    ╎    ╎   22     …sentials.jl:10; length
     ╎    ╎    ╎    ╎    ╎    ╎   45     …perators.jl:276; !=
   45╎    ╎    ╎    ╎    ╎    ╎    45     …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎  2      …c/matmul.jl:406; gemv!(y::Vector{Float…
     ╎    ╎    ╎    ╎    ╎    ╎   2      …perators.jl:276; !=
    2╎    ╎    ╎    ╎    ╎    ╎    2      …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎  8      …c/matmul.jl:409; gemv!(y::Vector{Float…
    8╎    ╎    ╎    ╎    ╎    ╎   8      …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎  14     …c/matmul.jl:410; gemv!(y::Vector{Float…
     ╎    ╎    ╎    ╎    ╎    ╎   14     …romotion.jl:399; promote
     ╎    ╎    ╎    ╎    ╎    ╎    14     …romotion.jl:377; _promote
     ╎    ╎    ╎    ╎    ╎    ╎     14     …e/number.jl:7; convert
   14╎    ╎    ╎    ╎    ╎    ╎    ╎ 14     …se/float.jl:165; Float64
     ╎    ╎    ╎    ╎    ╎    ╎  14     …c/matmul.jl:411; gemv!(y::Vector{Float…
   14╎    ╎    ╎    ╎    ╎    ╎   14     …perators.jl:1296; in
   51╎    ╎    ╎    ╎    ╎    ╎  1739   …c/matmul.jl:415; gemv!(y::Vector{Float…
   64╎    ╎    ╎    ╎    ╎    ╎   64     …src/blas.jl:643; gemv!(trans::Char, a…
     ╎    ╎    ╎    ╎    ╎    ╎   15     …src/blas.jl:647; gemv!(trans::Char, a…
   15╎    ╎    ╎    ╎    ╎    ╎    15     …se/array.jl:190; size
   28╎    ╎    ╎    ╎    ╎    ╎   30     …src/blas.jl:648; gemv!(trans::Char, a…
    1╎    ╎    ╎    ╎    ╎    ╎    1      …sentials.jl:10; length
     ╎    ╎    ╎    ╎    ╎    ╎    1      …perators.jl:276; !=
    1╎    ╎    ╎    ╎    ╎    ╎     1      …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎   2      …src/blas.jl:657; gemv!(trans::Char, a…
     ╎    ╎    ╎    ╎    ╎    ╎    2      …src/blas.jl:177; vec_pointer_stride
     ╎    ╎    ╎    ╎    ╎    ╎     2      …actarray.jl:1237; pointer
    2╎    ╎    ╎    ╎    ╎    ╎    ╎ 2      …/pointer.jl:65; unsafe_convert
     ╎    ╎    ╎    ╎    ╎    ╎   16     …src/blas.jl:659; gemv!(trans::Char, a…
     ╎    ╎    ╎    ╎    ╎    ╎    16     …actarray.jl:1237; pointer
   16╎    ╎    ╎    ╎    ╎    ╎     16     …/pointer.jl:65; unsafe_convert
     ╎    ╎    ╎    ╎    ╎    ╎   22     …src/blas.jl:666; gemv!(trans::Char, a…
     ╎    ╎    ╎    ╎    ╎    ╎    22     …perators.jl:587; max
     ╎    ╎    ╎    ╎    ╎    ╎     22     …romotion.jl:532; max
   12╎    ╎    ╎    ╎    ╎    ╎    ╎ 12     …sentials.jl:647; ifelse
   10╎    ╎    ╎    ╎    ╎    ╎    ╎ 10     @Base/int.jl:83; <
 1504╎    ╎    ╎    ╎    ╎    ╎   1539   …src/blas.jl:667; gemv!(trans::Char, a…
     ╎    ╎    ╎    ╎    ╎    ╎    35     …sentials.jl:543; cconvert
     ╎    ╎    ╎    ╎    ╎    ╎     35     …fpointer.jl:105; convert
   32╎    ╎    ╎    ╎    ╎    ╎    ╎ 35     …refvalue.jl:8; RefValue
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  3      …ase/char.jl:185; convert
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   3      …ase/char.jl:175; Int8
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    3      @Base/int.jl:538; >>>
    3╎    ╎    ╎    ╎    ╎    ╎    ╎     3      @Base/int.jl:530; >>>
   67╎    ╎    ╎    ╎    ╎  2869   …pse_comp.jl:75; eclipse_compute_quantities!…
   11╎    ╎    ╎    ╎    ╎   11     …se/array.jl:0; setindex!
   39╎    ╎    ╎    ╎    ╎   39     …se/array.jl:1021; setindex!
    4╎    ╎    ╎    ╎    ╎   4      …se/array.jl:0; vcat(::Vector{Float64}, ::V…
   59╎    ╎    ╎    ╎    ╎   59     …se/array.jl:2031; vcat(::Vector{Float64}, …
     ╎    ╎    ╎    ╎    ╎   12     …se/array.jl:2033; vcat(::Vector{Float64}, …
     ╎    ╎    ╎    ╎    ╎    12     …se/tuple.jl:72; iterate
     ╎    ╎    ╎    ╎    ╎     12     …se/tuple.jl:72; iterate
   12╎    ╎    ╎    ╎    ╎    ╎ 12     …se/tuple.jl:31; getindex
     ╎    ╎    ╎    ╎    ╎   4      …se/array.jl:2034; vcat(::Vector{Float64}, …
    1╎    ╎    ╎    ╎    ╎    1      …sentials.jl:10; length
    3╎    ╎    ╎    ╎    ╎    3      @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎   10     …se/array.jl:2035; vcat(::Vector{Float64}, …
     ╎    ╎    ╎    ╎    ╎    10     …se/tuple.jl:72; iterate
   10╎    ╎    ╎    ╎    ╎     10     …se/tuple.jl:31; getindex
     ╎    ╎    ╎    ╎    ╎   2162   …se/array.jl:2036; vcat(::Vector{Float64}, …
 2162╎    ╎    ╎    ╎    ╎    2162   …ase/boot.jl:477; Array
     ╎    ╎    ╎    ╎    ╎   14     …se/array.jl:2039; vcat(::Vector{Float64}, …
   14╎    ╎    ╎    ╎    ╎    14     …sentials.jl:10; length
   16╎    ╎    ╎    ╎    ╎   50     …se/array.jl:2040; vcat(::Vector{Float64}, …
    5╎    ╎    ╎    ╎    ╎    5      …sentials.jl:10; length
   29╎    ╎    ╎    ╎    ╎    29     @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎   390    …se/array.jl:2041; vcat(::Vector{Float64}, …
     ╎    ╎    ╎    ╎    ╎    20     …se/array.jl:331; unsafe_copyto!
     ╎    ╎    ╎    ╎    ╎     20     …actarray.jl:1240; pointer
    3╎    ╎    ╎    ╎    ╎    ╎ 3      …/pointer.jl:282; +
   17╎    ╎    ╎    ╎    ╎    ╎ 17     …/pointer.jl:65; unsafe_convert
     ╎    ╎    ╎    ╎    ╎    13     …se/array.jl:332; unsafe_copyto!
     ╎    ╎    ╎    ╎    ╎     13     …actarray.jl:1240; pointer
   13╎    ╎    ╎    ╎    ╎    ╎ 13     …/pointer.jl:65; unsafe_convert
     ╎    ╎    ╎    ╎    ╎    357    …se/array.jl:337; unsafe_copyto!
  336╎    ╎    ╎    ╎    ╎     347    …ase/cmem.jl:26; memmove
     ╎    ╎    ╎    ╎    ╎    ╎ 11     …sentials.jl:543; cconvert
     ╎    ╎    ╎    ╎    ╎    ╎  11     …e/number.jl:7; convert
     ╎    ╎    ╎    ╎    ╎    ╎   11     …ase/boot.jl:789; UInt64
     ╎    ╎    ╎    ╎    ╎    ╎    11     …ase/boot.jl:759; toUInt64
     ╎    ╎    ╎    ╎    ╎    ╎     11     …ase/boot.jl:648; check_top_bit
   11╎    ╎    ╎    ╎    ╎    ╎    ╎ 11     …ase/boot.jl:638; is_top_bit_set
   10╎    ╎    ╎    ╎    ╎     10     @Base/int.jl:88; *
     ╎    ╎    ╎    ╎    ╎   1      …se/array.jl:2043; vcat(::Vector{Float64}, …
    1╎    ╎    ╎    ╎    ╎    1      …se/tuple.jl:72; iterate
   17╎    ╎    ╎    ╎    ╎   17     …se/array.jl:2044; vcat(::Vector{Float64}, …
    4╎    ╎    ╎    ╎    ╎   4      …sentials.jl:0; getindex
   25╎    ╎    ╎    ╎    ╎   25     …sentials.jl:13; getindex
    1╎    ╎    ╎    ╎    ╎  6      …pse_comp.jl:76; eclipse_compute_quantities!…
    5╎    ╎    ╎    ╎    ╎   5      …se/range.jl:901; iterate
     ╎    ╎    ╎    ╎    ╎  2937   …pse_comp.jl:80; eclipse_compute_quantities!…
    2╎    ╎    ╎    ╎    ╎   2      …se/array.jl:0; setindex!
   51╎    ╎    ╎    ╎    ╎   51     …se/array.jl:1021; setindex!
     ╎    ╎    ╎    ╎    ╎   2881   …roadcast.jl:903; materialize
     ╎    ╎    ╎    ╎    ╎    2820   …roadcast.jl:928; copy
     ╎    ╎    ╎    ╎    ╎     491    …roadcast.jl:956; copyto!
    5╎    ╎    ╎    ╎    ╎    ╎ 5      …roadcast.jl:0; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎ 63     …roadcast.jl:1000; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎  63     …roadcast.jl:983; preprocess
     ╎    ╎    ╎    ╎    ╎    ╎   63     …roadcast.jl:986; preprocess_args
   22╎    ╎    ╎    ╎    ╎    ╎    22     …se/array.jl:0; size
     ╎    ╎    ╎    ╎    ╎    ╎    19     …roadcast.jl:984; preprocess
     ╎    ╎    ╎    ╎    ╎    ╎     19     …roadcast.jl:977; broadcast_unalias
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 19     …actarray.jl:1481; unalias
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  19     …actarray.jl:1516; mightalias
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   18     …actarray.jl:1523; _isdisjoint
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    18     …perators.jl:276; !=
   18╎    ╎    ╎    ╎    ╎    ╎    ╎     18     …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …actarray.jl:1540; dataids
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …actarray.jl:1237; pointer
    1╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …/pointer.jl:65; unsafe_convert
     ╎    ╎    ╎    ╎    ╎    ╎    22     …roadcast.jl:987; preprocess_args
     ╎    ╎    ╎    ╎    ╎    ╎     22     …roadcast.jl:984; preprocess
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 22     …roadcast.jl:977; broadcast_unalias
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  22     …actarray.jl:1481; unalias
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   22     …actarray.jl:1516; mightalias
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    22     …actarray.jl:1523; _isdisjoint
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     22     …perators.jl:276; !=
   22╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 22     …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎ 423    …roadcast.jl:1003; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎  1      …simdloop.jl:72; macro expansion
    1╎    ╎    ╎    ╎    ╎    ╎   1      @Base/int.jl:83; <
   47╎    ╎    ╎    ╎    ╎    ╎  47     …simdloop.jl:75; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎  345    …simdloop.jl:77; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎   345    …roadcast.jl:1004; macro expansion
  296╎    ╎    ╎    ╎    ╎    ╎    296    …se/array.jl:1021; setindex!
     ╎    ╎    ╎    ╎    ╎    ╎    49     …roadcast.jl:636; getindex
     ╎    ╎    ╎    ╎    ╎    ╎     47     …roadcast.jl:681; _broadcast_getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 47     …roadcast.jl:705; _getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  47     …roadcast.jl:675; _broadcast_getin…
   47╎    ╎    ╎    ╎    ╎    ╎    ╎   47     …sentials.jl:13; getindex
     ╎    ╎    ╎    ╎    ╎    ╎     2      …roadcast.jl:682; _broadcast_getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 2      …roadcast.jl:709; _broadcast_getind…
    2╎    ╎    ╎    ╎    ╎    ╎    ╎  2      …se/float.jl:409; +
     ╎    ╎    ╎    ╎    ╎    ╎  30     …simdloop.jl:78; macro expansion
   30╎    ╎    ╎    ╎    ╎    ╎   30     @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎     2329   …roadcast.jl:223; similar
     ╎    ╎    ╎    ╎    ╎    ╎ 2329   …roadcast.jl:224; similar
     ╎    ╎    ╎    ╎    ╎    ╎  2329   …actarray.jl:876; similar
     ╎    ╎    ╎    ╎    ╎    ╎   2329   …actarray.jl:877; similar
     ╎    ╎    ╎    ╎    ╎    ╎    2329   …ase/boot.jl:494; Array
     ╎    ╎    ╎    ╎    ╎    ╎     2329   …ase/boot.jl:486; Array
 2329╎    ╎    ╎    ╎    ╎    ╎    ╎ 2329   …ase/boot.jl:477; Array
     ╎    ╎    ╎    ╎    ╎    61     …roadcast.jl:306; instantiate
     ╎    ╎    ╎    ╎    ╎     61     …roadcast.jl:524; combine_axes
     ╎    ╎    ╎    ╎    ╎    ╎ 23     …actarray.jl:98; axes
   23╎    ╎    ╎    ╎    ╎    ╎  23     …se/array.jl:191; size
     ╎    ╎    ╎    ╎    ╎    ╎ 38     …roadcast.jl:543; broadcast_shape
   16╎    ╎    ╎    ╎    ╎    ╎  16     …roadcast.jl:0; _bcs
    8╎    ╎    ╎    ╎    ╎    ╎  22     …roadcast.jl:549; _bcs
     ╎    ╎    ╎    ╎    ╎    ╎   14     …roadcast.jl:555; _bcs1
   14╎    ╎    ╎    ╎    ╎    ╎    14     …roadcast.jl:557; _bcsm
    3╎    ╎    ╎    ╎    ╎   3      …sentials.jl:13; getindex
    2╎    ╎    ╎    ╎    ╎  5      …pse_comp.jl:81; eclipse_compute_quantities!…
    3╎    ╎    ╎    ╎    ╎   3      …se/range.jl:901; iterate
     ╎    ╎    ╎    ╎    ╎  2096   …pse_comp.jl:84; eclipse_compute_quantities!…
    9╎    ╎    ╎    ╎    ╎   9      …geometry.jl:96; calc_mu(xyz::SubArray{Floa…
   80╎    ╎    ╎    ╎    ╎   2086   …geometry.jl:108; calc_mu_grid!(A::Matrix{V…
    1╎    ╎    ╎    ╎    ╎    1      …se/array.jl:1021; setindex!
   11╎    ╎    ╎    ╎    ╎    11     …sentials.jl:0; getindex
     ╎    ╎    ╎    ╎    ╎    23     …subarray.jl:184; view
   14╎    ╎    ╎    ╎    ╎     23     …actarray.jl:702; checkbounds
     ╎    ╎    ╎    ╎    ╎    ╎ 9      …actarray.jl:687; checkbounds
     ╎    ╎    ╎    ╎    ╎    ╎  9      …actarray.jl:768; checkindex
     ╎    ╎    ╎    ╎    ╎    ╎   9      …actarray.jl:763; checkindex
    9╎    ╎    ╎    ╎    ╎    ╎    9      @Base/int.jl:513; <
   23╎    ╎    ╎    ╎    ╎    23     …/generic.jl:595; norm(itr::SubArray{Float…
   52╎    ╎    ╎    ╎    ╎    52     …geometry.jl:95; calc_mu(xyz::SubArray{Flo…
   17╎    ╎    ╎    ╎    ╎    1896   …geometry.jl:96; calc_mu(xyz::SubArray{Flo…
    9╎    ╎    ╎    ╎    ╎     9      …se/float.jl:411; *
  168╎    ╎    ╎    ╎    ╎     168    …se/float.jl:412; /
   32╎    ╎    ╎    ╎    ╎     1436   …/generic.jl:596; norm
   74╎    ╎    ╎    ╎    ╎    ╎ 74     …educedim.jl:0; norm(itr::SubArray{Float…
   16╎    ╎    ╎    ╎    ╎    ╎ 16     …/generic.jl:0; generic_norm2(x::SubArra…
   24╎    ╎    ╎    ╎    ╎    ╎ 24     …/generic.jl:462; generic_norm2(x::SubAr…
  111╎    ╎    ╎    ╎    ╎    ╎ 111    …/generic.jl:595; norm(itr::SubArray{Flo…
     ╎    ╎    ╎    ╎    ╎    ╎ 54     …/generic.jl:596; norm(itr::SubArray{Flo…
     ╎    ╎    ╎    ╎    ╎    ╎  54     …actarray.jl:1220; isempty
     ╎    ╎    ╎    ╎    ╎    ╎   40     …actarray.jl:315; length
     ╎    ╎    ╎    ╎    ╎    ╎    40     …subarray.jl:63; size
     ╎    ╎    ╎    ╎    ╎    ╎     40     …subarray.jl:490; axes
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 40     …subarray.jl:495; _indices_sub
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  40     …se/range.jl:706; axes
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   40     …se/range.jl:761; length
   40╎    ╎    ╎    ╎    ╎    ╎    ╎    40     @Base/int.jl:86; -
   14╎    ╎    ╎    ╎    ╎    ╎   14     …romotion.jl:521; ==
   11╎    ╎    ╎    ╎    ╎    ╎ 11     …/generic.jl:597; norm(itr::SubArray{Flo…
     ╎    ╎    ╎    ╎    ╎    ╎ 1114   …/generic.jl:598; norm(itr::SubArray{Flo…
    2╎    ╎    ╎    ╎    ╎    ╎  1114   …rc/dense.jl:106; norm2
     ╎    ╎    ╎    ╎    ╎    ╎   27     …actarray.jl:315; length
     ╎    ╎    ╎    ╎    ╎    ╎    27     …subarray.jl:63; size
     ╎    ╎    ╎    ╎    ╎    ╎     27     …subarray.jl:490; axes
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 27     …subarray.jl:495; _indices_sub
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  27     …se/range.jl:706; axes
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   12     …se/range.jl:761; length
   12╎    ╎    ╎    ╎    ╎    ╎    ╎    12     @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   15     …se/range.jl:469; oneto
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    15     …se/range.jl:467; OneTo
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     15     …se/range.jl:454; OneTo
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 15     …romotion.jl:532; max
   15╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  15     …sentials.jl:647; ifelse
   15╎    ╎    ╎    ╎    ╎    ╎   15     @Base/int.jl:83; <
   22╎    ╎    ╎    ╎    ╎    ╎   22     …e/reduce.jl:0; _mapreduce(f::typeof(L…
   21╎    ╎    ╎    ╎    ╎    ╎   21     …e/reduce.jl:428; _mapreduce(f::typeof…
   96╎    ╎    ╎    ╎    ╎    ╎   96     …/generic.jl:0; generic_norm2(x::SubAr…
   64╎    ╎    ╎    ╎    ╎    ╎   64     …/generic.jl:462; generic_norm2(x::Sub…
     ╎    ╎    ╎    ╎    ╎    ╎   530    …/generic.jl:463; generic_norm2(x::Sub…
     ╎    ╎    ╎    ╎    ╎    ╎    530    …/generic.jl:527; normInf
     ╎    ╎    ╎    ╎    ╎    ╎     530    …/generic.jl:453; generic_normInf
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 530    …educedim.jl:357; mapreduce
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  530    …educedim.jl:357; #mapreduce#821
    1╎    ╎    ╎    ╎    ╎    ╎    ╎   530    …educedim.jl:365; _mapreduce_dim
   81╎    ╎    ╎    ╎    ╎    ╎    ╎    81     …e/reduce.jl:0; _mapreduce(f::ty…
  133╎    ╎    ╎    ╎    ╎    ╎    ╎    133    …e/reduce.jl:428; _mapreduce(f::…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    10     …e/reduce.jl:429; _mapreduce(f::…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     10     …/indices.jl:486; LinearIndices
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 10     …subarray.jl:490; axes
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  10     …subarray.jl:495; _indices_sub
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   10     …se/range.jl:706; axes
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    10     …se/range.jl:761; length
   10╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     10     @Base/int.jl:86; -
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    11     …e/reduce.jl:431; _mapreduce(f::…
   11╎    ╎    ╎    ╎    ╎    ╎    ╎     11     …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    24     …e/reduce.jl:438; _mapreduce(f::…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     24     …subarray.jl:323; getindex
   24╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 24     …ase/Base.jl:37; getproperty
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    79     …e/reduce.jl:440; _mapreduce(f::…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     2      …ase/math.jl:908; max
    2╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 2      …se/float.jl:410; -
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     34     …ase/math.jl:909; max
   27╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 27     …sentials.jl:647; ifelse
    7╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 7      …oatfuncs.jl:15; signbit
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     20     …ase/math.jl:911; max
   20╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 20     …sentials.jl:647; ifelse
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     23     …/generic.jl:639; norm
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 23     …/generic.jl:639; norm
   23╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  23     …se/float.jl:610; abs
  111╎    ╎    ╎    ╎    ╎    ╎    ╎    114    …e/reduce.jl:441; _mapreduce(f::…
    3╎    ╎    ╎    ╎    ╎    ╎    ╎     3      @Base/int.jl:83; <
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    77     …e/reduce.jl:443; _mapreduce(f::…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     26     …ase/math.jl:909; max
   26╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 26     …sentials.jl:647; ifelse
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     31     …ase/math.jl:911; max
   31╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 31     …sentials.jl:647; ifelse
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     20     …/generic.jl:639; norm
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 20     …/generic.jl:639; norm
   20╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  20     …se/float.jl:610; abs
     ╎    ╎    ╎    ╎    ╎    ╎   56     …/generic.jl:464; generic_norm2(x::Sub…
   31╎    ╎    ╎    ╎    ╎    ╎    56     …se/float.jl:635; isinf
   25╎    ╎    ╎    ╎    ╎    ╎     25     …se/float.jl:610; abs
     ╎    ╎    ╎    ╎    ╎    ╎   80     …/generic.jl:465; generic_norm2(x::Sub…
     ╎    ╎    ╎    ╎    ╎    ╎    80     …actarray.jl:1215; iterate
     ╎    ╎    ╎    ╎    ╎    ╎     3      …actarray.jl:321; eachindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 3      …actarray.jl:137; axes1
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  3      …subarray.jl:490; axes
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   3      …subarray.jl:495; _indices_sub
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    3      …se/range.jl:706; axes
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     3      …se/range.jl:761; length
    3╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 3      @Base/int.jl:86; -
   29╎    ╎    ╎    ╎    ╎    ╎     29     …actarray.jl:0; ifelse
     ╎    ╎    ╎    ╎    ╎    ╎     24     …actarray.jl:1215; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 24     …se/range.jl:897; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  24     …se/range.jl:672; isempty
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   24     …perators.jl:378; >
   24╎    ╎    ╎    ╎    ╎    ╎    ╎    24     @Base/int.jl:83; <
     ╎    ╎    ╎    ╎    ╎    ╎     24     …actarray.jl:1217; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 24     …subarray.jl:323; getindex
    1╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …ase/Base.jl:37; getproperty
   23╎    ╎    ╎    ╎    ╎    ╎    ╎  23     …sentials.jl:13; getindex
   21╎    ╎    ╎    ╎    ╎    ╎   60     …/generic.jl:467; generic_norm2(x::Sub…
     ╎    ╎    ╎    ╎    ╎    ╎    11     …se/float.jl:623; isfinite
    3╎    ╎    ╎    ╎    ╎    ╎     3      …se/float.jl:410; -
     ╎    ╎    ╎    ╎    ╎    ╎     8      …se/float.jl:620; isnan
    8╎    ╎    ╎    ╎    ╎    ╎    ╎ 8      …se/float.jl:535; !=
     ╎    ╎    ╎    ╎    ╎    ╎    28     …perators.jl:587; *
   28╎    ╎    ╎    ╎    ╎    ╎     28     …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎   27     …/generic.jl:468; generic_norm2(x::Sub…
     ╎    ╎    ╎    ╎    ╎    ╎    27     …/generic.jl:459; norm_sqr
     ╎    ╎    ╎    ╎    ╎    ╎     27     …e/number.jl:189; abs2
   27╎    ╎    ╎    ╎    ╎    ╎    ╎ 27     …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎   73     …/generic.jl:470; generic_norm2(x::Sub…
     ╎    ╎    ╎    ╎    ╎    ╎    25     …actarray.jl:1215; iterate
     ╎    ╎    ╎    ╎    ╎    ╎     25     …se/range.jl:901; iterate
   25╎    ╎    ╎    ╎    ╎    ╎    ╎ 25     …romotion.jl:521; ==
   20╎    ╎    ╎    ╎    ╎    ╎    20     …actarray.jl:1216; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    28     …actarray.jl:1217; iterate
     ╎    ╎    ╎    ╎    ╎    ╎     28     …subarray.jl:323; getindex
   28╎    ╎    ╎    ╎    ╎    ╎    ╎ 28     …sentials.jl:13; getindex
     ╎    ╎    ╎    ╎    ╎    ╎   22     …/generic.jl:473; generic_norm2(x::Sub…
   22╎    ╎    ╎    ╎    ╎    ╎    22     …se/float.jl:409; +
     ╎    ╎    ╎    ╎    ╎    ╎   19     …/generic.jl:476; generic_norm2(x::Sub…
   19╎    ╎    ╎    ╎    ╎    ╎    19     …ase/math.jl:686; sqrt
     ╎    ╎    ╎    ╎    ╎     266    …c/matmul.jl:15; dot
    1╎    ╎    ╎    ╎    ╎    ╎ 1      …src/blas.jl:0; dot
     ╎    ╎    ╎    ╎    ╎    ╎ 54     …src/blas.jl:393; dot
     ╎    ╎    ╎    ╎    ╎    ╎  54     …actarray.jl:315; length
     ╎    ╎    ╎    ╎    ╎    ╎   54     …subarray.jl:63; size
     ╎    ╎    ╎    ╎    ╎    ╎    54     …subarray.jl:490; axes
     ╎    ╎    ╎    ╎    ╎    ╎     54     …subarray.jl:495; _indices_sub
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 54     …se/range.jl:706; axes
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  23     …se/range.jl:761; length
    1╎    ╎    ╎    ╎    ╎    ╎    ╎   1      @Base/int.jl:87; +
   22╎    ╎    ╎    ╎    ╎    ╎    ╎   22     @Base/int.jl:86; -
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  31     …se/range.jl:469; oneto
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   31     …se/range.jl:467; OneTo
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    31     …se/range.jl:454; OneTo
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     31     …romotion.jl:532; max
   31╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 31     …sentials.jl:647; ifelse
   11╎    ╎    ╎    ╎    ╎    ╎ 211    …src/blas.jl:395; dot
  175╎    ╎    ╎    ╎    ╎    ╎  175    …src/blas.jl:345; dot
     ╎    ╎    ╎    ╎    ╎    ╎  25     …src/blas.jl:177; vec_pointer_stride
     ╎    ╎    ╎    ╎    ╎    ╎   25     …src/blas.jl:177; vec_pointer_stride
     ╎    ╎    ╎    ╎    ╎    ╎    25     …actarray.jl:1237; pointer
     ╎    ╎    ╎    ╎    ╎    ╎     25     …subarray.jl:472; unsafe_convert
   25╎    ╎    ╎    ╎    ╎    ╎    ╎ 25     …/pointer.jl:282; +
    1╎    ╎    ╎    ╎    ╎   1      …geometry.jl:109; calc_mu_grid!(A::Matrix{V…
     ╎    ╎    ╎    ╎    ╎  163    …pse_comp.jl:86; eclipse_compute_quantities!…
     ╎    ╎    ╎    ╎    ╎   163    …roadcast.jl:903; materialize
     ╎    ╎    ╎    ╎    ╎    161    …roadcast.jl:928; copy
     ╎    ╎    ╎    ╎    ╎     143    …roadcast.jl:956; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎ 90     …roadcast.jl:1015; copyto!
   90╎    ╎    ╎    ╎    ╎    ╎  90     …ase/boot.jl:477; Array
     ╎    ╎    ╎    ╎    ╎    ╎ 6      …roadcast.jl:1018; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎  6      …roadcast.jl:983; preprocess
     ╎    ╎    ╎    ╎    ╎    ╎   6      …roadcast.jl:986; preprocess_args
     ╎    ╎    ╎    ╎    ╎    ╎    6      …roadcast.jl:984; preprocess
     ╎    ╎    ╎    ╎    ╎    ╎     6      …roadcast.jl:977; broadcast_unalias
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 6      …actarray.jl:1481; unalias
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  6      …actarray.jl:1516; mightalias
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   6      …actarray.jl:1539; dataids
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    6      …flection.jl:611; objectid
    6╎    ╎    ╎    ╎    ╎    ╎    ╎     6      …flection.jl:617; _objectid
     ╎    ╎    ╎    ╎    ╎    ╎ 2      …roadcast.jl:1019; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎  2      …terators.jl:1271; partition
     ╎    ╎    ╎    ╎    ╎    ╎   2      …terators.jl:1279; PartitionIterator
     ╎    ╎    ╎    ╎    ╎    ╎    2      …rraymath.jl:41; vec
     ╎    ╎    ╎    ╎    ╎    ╎     2      …pedarray.jl:117; reshape
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 2      …pedarray.jl:112; reshape
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  2      …pedarray.jl:178; _reshape
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   2      …pedarray.jl:193; __reshape
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    2      …se/tuple.jl:291; map
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     2      …inverses.jl:90; SignedMultipli…
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …inverses.jl:82; Base.Multipli…
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …inverses.jl:87; Base.Multipli…
     ╎    ╎    ╎    ╎    ╎    ╎ 29     …roadcast.jl:1021; copyto!
    2╎    ╎    ╎    ╎    ╎    ╎  2      …simdloop.jl:0; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎  1      …simdloop.jl:70; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎   1      …ensional.jl:661; simd_outer_range
     ╎    ╎    ╎    ╎    ╎    ╎    1      …pedarray.jl:220; ind2sub_rs
     ╎    ╎    ╎    ╎    ╎    ╎     1      …pedarray.jl:223; _ind2sub_rs
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …inverses.jl:172; divrem
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …inverses.jl:160; div
    1╎    ╎    ╎    ╎    ╎    ╎    ╎   1      @Base/int.jl:88; *
   10╎    ╎    ╎    ╎    ╎    ╎  11     …simdloop.jl:75; macro expansion
    1╎    ╎    ╎    ╎    ╎    ╎   1      @Base/int.jl:83; <
     ╎    ╎    ╎    ╎    ╎    ╎  12     …simdloop.jl:77; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎   12     …roadcast.jl:1022; macro expansion
    5╎    ╎    ╎    ╎    ╎    ╎    5      …se/array.jl:1021; setindex!
     ╎    ╎    ╎    ╎    ╎    ╎    7      …roadcast.jl:636; getindex
     ╎    ╎    ╎    ╎    ╎    ╎     7      …roadcast.jl:682; _broadcast_getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 7      …roadcast.jl:709; _broadcast_getind…
    7╎    ╎    ╎    ╎    ╎    ╎    ╎  7      …se/float.jl:536; <
     ╎    ╎    ╎    ╎    ╎    ╎  3      …simdloop.jl:84; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎   3      …enerator.jl:47; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    3      none:0; #30
    1╎    ╎    ╎    ╎    ╎    ╎     1      …ensional.jl:665; skip_len_I
     ╎    ╎    ╎    ╎    ╎    ╎     1      …ensional.jl:666; skip_len_I
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …ensional.jl:108; ==
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …se/tuple.jl:482; ==
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …se/tuple.jl:486; _eq
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎     1      …ensional.jl:667; skip_len_I
    1╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎    ╎ 1      …roadcast.jl:1025; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎  1      …simdloop.jl:77; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎   1      …roadcast.jl:1026; macro expansion
    1╎    ╎    ╎    ╎    ╎    ╎    1      …se/array.jl:1021; setindex!
     ╎    ╎    ╎    ╎    ╎    ╎ 14     …roadcast.jl:1028; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎  14     …bitarray.jl:356; dumpbitcache
    2╎    ╎    ╎    ╎    ╎    ╎   2      …bitarray.jl:0; copy_to_bitarray_chunk…
    2╎    ╎    ╎    ╎    ╎    ╎   2      …bitarray.jl:283; copy_to_bitarray_chu…
     ╎    ╎    ╎    ╎    ╎    ╎   2      …bitarray.jl:284; copy_to_bitarray_chu…
     ╎    ╎    ╎    ╎    ╎    ╎    2      …bitarray.jl:126; get_chunks_id
    2╎    ╎    ╎    ╎    ╎    ╎     2      @Base/int.jl:86; -
     ╎    ╎    ╎    ╎    ╎    ╎   4      …bitarray.jl:320; copy_to_bitarray_chu…
     ╎    ╎    ╎    ╎    ╎    ╎    2      …bitarray.jl:276; pack8bools
     ╎    ╎    ╎    ╎    ╎    ╎     1      @Base/int.jl:538; >>>
    1╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      @Base/int.jl:530; >>>
    1╎    ╎    ╎    ╎    ╎    ╎     1      @Base/int.jl:372; |
     ╎    ╎    ╎    ╎    ╎    ╎    1      …bitarray.jl:277; pack8bools
     ╎    ╎    ╎    ╎    ╎    ╎     1      @Base/int.jl:538; >>>
    1╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      @Base/int.jl:530; >>>
     ╎    ╎    ╎    ╎    ╎    ╎    1      …bitarray.jl:278; pack8bools
     ╎    ╎    ╎    ╎    ╎    ╎     1      @Base/int.jl:538; >>>
    1╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      @Base/int.jl:530; >>>
     ╎    ╎    ╎    ╎    ╎    ╎   2      …bitarray.jl:325; copy_to_bitarray_chu…
    2╎    ╎    ╎    ╎    ╎    ╎    2      …se/range.jl:901; iterate
    1╎    ╎    ╎    ╎    ╎    ╎   1      …bitarray.jl:344; copy_to_bitarray_chu…
     ╎    ╎    ╎    ╎    ╎    ╎   1      @Base/int.jl:536; <<
    1╎    ╎    ╎    ╎    ╎    ╎    1      @Base/int.jl:529; <<
    1╎    ╎    ╎    ╎    ╎    ╎ 1      @Base/int.jl:0; copyto!
     ╎    ╎    ╎    ╎    ╎     18     …roadcast.jl:223; similar
     ╎    ╎    ╎    ╎    ╎    ╎ 18     …roadcast.jl:226; similar
     ╎    ╎    ╎    ╎    ╎    ╎  18     …actarray.jl:876; similar
     ╎    ╎    ╎    ╎    ╎    ╎   18     …actarray.jl:877; similar
     ╎    ╎    ╎    ╎    ╎    ╎    18     …bitarray.jl:71; BitArray
    1╎    ╎    ╎    ╎    ╎    ╎     1      …bitarray.jl:28; BitMatrix(::UndefIn…
     ╎    ╎    ╎    ╎    ╎    ╎     14     …bitarray.jl:37; BitMatrix(::UndefIn…
   14╎    ╎    ╎    ╎    ╎    ╎    ╎ 14     …ase/boot.jl:477; Array
    3╎    ╎    ╎    ╎    ╎    ╎     3      …bitarray.jl:39; BitMatrix(::UndefIn…
     ╎    ╎    ╎    ╎    ╎    2      …roadcast.jl:306; instantiate
     ╎    ╎    ╎    ╎    ╎     2      …roadcast.jl:524; combine_axes
     ╎    ╎    ╎    ╎    ╎    ╎ 2      …roadcast.jl:543; broadcast_shape
    2╎    ╎    ╎    ╎    ╎    ╎  2      …roadcast.jl:547; _bcs
   32╎    ╎    ╎    ╎    ╎  4300   …pse_comp.jl:89; eclipse_compute_quantities!…
   21╎    ╎    ╎    ╎    ╎   21     …se/array.jl:371; _copyto_impl!(dest::Vecto…
    4╎    ╎    ╎    ╎    ╎   4      …sentials.jl:0; projected!(A::Matrix{Vector…
   14╎    ╎    ╎    ╎    ╎   14     …/generic.jl:595; norm(itr::Vector{Float64}…
     ╎    ╎    ╎    ╎    ╎   1126   …_physics.jl:69; projected!(A::Matrix{Vecto…
    2╎    ╎    ╎    ╎    ╎    2      …actarray.jl:0; getindex
     ╎    ╎    ╎    ╎    ╎    13     …se/array.jl:973; getindex
    1╎    ╎    ╎    ╎    ╎     1      …actarray.jl:700; checkbounds
    6╎    ╎    ╎    ╎    ╎     12     …actarray.jl:702; checkbounds
     ╎    ╎    ╎    ╎    ╎    ╎ 6      …actarray.jl:687; checkbounds
     ╎    ╎    ╎    ╎    ╎    ╎  6      …actarray.jl:768; checkindex
     ╎    ╎    ╎    ╎    ╎    ╎   6      …actarray.jl:763; checkindex
    6╎    ╎    ╎    ╎    ╎    ╎    6      @Base/int.jl:513; <
     ╎    ╎    ╎    ╎    ╎    945    …se/array.jl:975; getindex
     ╎    ╎    ╎    ╎    ╎     945    …actarray.jl:831; similar
     ╎    ╎    ╎    ╎    ╎    ╎ 945    …se/array.jl:420; similar
     ╎    ╎    ╎    ╎    ╎    ╎  945    …ase/boot.jl:486; Array
  945╎    ╎    ╎    ╎    ╎    ╎   945    …ase/boot.jl:477; Array
     ╎    ╎    ╎    ╎    ╎    162    …se/array.jl:977; getindex
   27╎    ╎    ╎    ╎    ╎     162    …se/array.jl:368; copyto!
   12╎    ╎    ╎    ╎    ╎    ╎ 12     …se/array.jl:0; _copyto_impl!(dest::Vect…
   14╎    ╎    ╎    ╎    ╎    ╎ 14     …se/array.jl:371; _copyto_impl!(dest::Ve…
     ╎    ╎    ╎    ╎    ╎    ╎ 1      …se/array.jl:372; _copyto_impl!(dest::Ve…
    1╎    ╎    ╎    ╎    ╎    ╎  1      …romotion.jl:521; ==
    8╎    ╎    ╎    ╎    ╎    ╎ 8      …se/array.jl:373; _copyto_impl!(dest::Ve…
     ╎    ╎    ╎    ╎    ╎    ╎ 18     …se/array.jl:374; _copyto_impl!(dest::Ve…
     ╎    ╎    ╎    ╎    ╎    ╎  14     …actarray.jl:702; checkbounds
     ╎    ╎    ╎    ╎    ╎    ╎   14     …actarray.jl:687; checkbounds
     ╎    ╎    ╎    ╎    ╎    ╎    4      …actarray.jl:768; checkindex
     ╎    ╎    ╎    ╎    ╎    ╎     4      …actarray.jl:763; checkindex
    4╎    ╎    ╎    ╎    ╎    ╎    ╎ 4      @Base/int.jl:513; <
     ╎    ╎    ╎    ╎    ╎    ╎    10     …actarray.jl:389; eachindex
     ╎    ╎    ╎    ╎    ╎    ╎     10     …actarray.jl:137; axes1
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 10     …actarray.jl:98; axes
   10╎    ╎    ╎    ╎    ╎    ╎    ╎  10     …se/array.jl:191; size
     ╎    ╎    ╎    ╎    ╎    ╎  4      …se/range.jl:5; Colon
     ╎    ╎    ╎    ╎    ╎    ╎   4      …se/range.jl:403; UnitRange
     ╎    ╎    ╎    ╎    ╎    ╎    4      …se/range.jl:414; unitrange_last
     ╎    ╎    ╎    ╎    ╎    ╎     4      …perators.jl:425; >=
    4╎    ╎    ╎    ╎    ╎    ╎    ╎ 4      @Base/int.jl:514; <=
     ╎    ╎    ╎    ╎    ╎    ╎ 31     …se/array.jl:375; _copyto_impl!(dest::Ve…
     ╎    ╎    ╎    ╎    ╎    ╎  16     …actarray.jl:702; checkbounds
     ╎    ╎    ╎    ╎    ╎    ╎   16     …actarray.jl:687; checkbounds
     ╎    ╎    ╎    ╎    ╎    ╎    16     …actarray.jl:768; checkindex
     ╎    ╎    ╎    ╎    ╎    ╎     16     …actarray.jl:763; checkindex
   16╎    ╎    ╎    ╎    ╎    ╎    ╎ 16     @Base/int.jl:513; <
    8╎    ╎    ╎    ╎    ╎    ╎  8      @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎    ╎  7      …se/range.jl:5; Colon
     ╎    ╎    ╎    ╎    ╎    ╎   7      …se/range.jl:403; UnitRange
    7╎    ╎    ╎    ╎    ╎    ╎    7      …se/range.jl:414; unitrange_last
     ╎    ╎    ╎    ╎    ╎    ╎ 51     …se/array.jl:376; _copyto_impl!(dest::Ve…
     ╎    ╎    ╎    ╎    ╎    ╎  6      …se/array.jl:332; unsafe_copyto!
     ╎    ╎    ╎    ╎    ╎    ╎   6      …actarray.jl:1240; pointer
    6╎    ╎    ╎    ╎    ╎    ╎    6      …/pointer.jl:282; +
     ╎    ╎    ╎    ╎    ╎    ╎  45     …se/array.jl:337; unsafe_copyto!
   39╎    ╎    ╎    ╎    ╎    ╎   45     …ase/cmem.jl:26; memmove
     ╎    ╎    ╎    ╎    ╎    ╎    6      …sentials.jl:543; cconvert
     ╎    ╎    ╎    ╎    ╎    ╎     6      …e/number.jl:7; convert
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 6      …ase/boot.jl:789; UInt64
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  6      …ase/boot.jl:759; toUInt64
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   6      …ase/boot.jl:648; check_top_bit
    6╎    ╎    ╎    ╎    ╎    ╎    ╎    6      …ase/boot.jl:638; is_top_bit_set
    4╎    ╎    ╎    ╎    ╎    4      …sentials.jl:13; getindex
     ╎    ╎    ╎    ╎    ╎   2900   …_physics.jl:70; projected!(A::Matrix{Vecto…
    2╎    ╎    ╎    ╎    ╎    2      …actarray.jl:0; getindex
     ╎    ╎    ╎    ╎    ╎    25     …se/array.jl:973; getindex
    4╎    ╎    ╎    ╎    ╎     4      …actarray.jl:700; checkbounds
   13╎    ╎    ╎    ╎    ╎     21     …actarray.jl:702; checkbounds
     ╎    ╎    ╎    ╎    ╎    ╎ 8      …actarray.jl:687; checkbounds
     ╎    ╎    ╎    ╎    ╎    ╎  8      …actarray.jl:768; checkindex
     ╎    ╎    ╎    ╎    ╎    ╎   8      …actarray.jl:763; checkindex
    8╎    ╎    ╎    ╎    ╎    ╎    8      @Base/int.jl:513; <
     ╎    ╎    ╎    ╎    ╎    1935   …se/array.jl:975; getindex
     ╎    ╎    ╎    ╎    ╎     1935   …actarray.jl:831; similar
     ╎    ╎    ╎    ╎    ╎    ╎ 1935   …se/array.jl:420; similar
     ╎    ╎    ╎    ╎    ╎    ╎  1935   …ase/boot.jl:486; Array
 1935╎    ╎    ╎    ╎    ╎    ╎   1935   …ase/boot.jl:477; Array
     ╎    ╎    ╎    ╎    ╎    242    …se/array.jl:977; getindex
   42╎    ╎    ╎    ╎    ╎     242    …se/array.jl:368; copyto!
    7╎    ╎    ╎    ╎    ╎    ╎ 7      …se/array.jl:0; _copyto_impl!(dest::Vect…
   18╎    ╎    ╎    ╎    ╎    ╎ 18     …se/array.jl:371; _copyto_impl!(dest::Ve…
     ╎    ╎    ╎    ╎    ╎    ╎ 4      …se/array.jl:372; _copyto_impl!(dest::Ve…
    4╎    ╎    ╎    ╎    ╎    ╎  4      …romotion.jl:521; ==
   20╎    ╎    ╎    ╎    ╎    ╎ 20     …se/array.jl:373; _copyto_impl!(dest::Ve…
     ╎    ╎    ╎    ╎    ╎    ╎ 39     …se/array.jl:374; _copyto_impl!(dest::Ve…
     ╎    ╎    ╎    ╎    ╎    ╎  33     …actarray.jl:702; checkbounds
     ╎    ╎    ╎    ╎    ╎    ╎   33     …actarray.jl:687; checkbounds
     ╎    ╎    ╎    ╎    ╎    ╎    18     …actarray.jl:768; checkindex
     ╎    ╎    ╎    ╎    ╎    ╎     18     …actarray.jl:763; checkindex
   18╎    ╎    ╎    ╎    ╎    ╎    ╎ 18     @Base/int.jl:513; <
     ╎    ╎    ╎    ╎    ╎    ╎    15     …actarray.jl:389; eachindex
     ╎    ╎    ╎    ╎    ╎    ╎     15     …actarray.jl:137; axes1
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 15     …actarray.jl:98; axes
   15╎    ╎    ╎    ╎    ╎    ╎    ╎  15     …se/array.jl:191; size
     ╎    ╎    ╎    ╎    ╎    ╎  6      …se/range.jl:5; Colon
     ╎    ╎    ╎    ╎    ╎    ╎   6      …se/range.jl:403; UnitRange
     ╎    ╎    ╎    ╎    ╎    ╎    6      …se/range.jl:414; unitrange_last
     ╎    ╎    ╎    ╎    ╎    ╎     6      …perators.jl:425; >=
    6╎    ╎    ╎    ╎    ╎    ╎    ╎ 6      @Base/int.jl:514; <=
     ╎    ╎    ╎    ╎    ╎    ╎ 38     …se/array.jl:375; _copyto_impl!(dest::Ve…
     ╎    ╎    ╎    ╎    ╎    ╎  15     …actarray.jl:702; checkbounds
     ╎    ╎    ╎    ╎    ╎    ╎   15     …actarray.jl:687; checkbounds
     ╎    ╎    ╎    ╎    ╎    ╎    15     …actarray.jl:768; checkindex
     ╎    ╎    ╎    ╎    ╎    ╎     15     …actarray.jl:763; checkindex
   15╎    ╎    ╎    ╎    ╎    ╎    ╎ 15     @Base/int.jl:513; <
   10╎    ╎    ╎    ╎    ╎    ╎  10     @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎    ╎  13     …se/range.jl:5; Colon
     ╎    ╎    ╎    ╎    ╎    ╎   13     …se/range.jl:403; UnitRange
   12╎    ╎    ╎    ╎    ╎    ╎    13     …se/range.jl:414; unitrange_last
     ╎    ╎    ╎    ╎    ╎    ╎     1      …perators.jl:425; >=
    1╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      @Base/int.jl:514; <=
     ╎    ╎    ╎    ╎    ╎    ╎ 74     …se/array.jl:376; _copyto_impl!(dest::Ve…
     ╎    ╎    ╎    ╎    ╎    ╎  1      …se/array.jl:331; unsafe_copyto!
     ╎    ╎    ╎    ╎    ╎    ╎   1      …actarray.jl:1240; pointer
    1╎    ╎    ╎    ╎    ╎    ╎    1      …/pointer.jl:282; +
     ╎    ╎    ╎    ╎    ╎    ╎  9      …se/array.jl:332; unsafe_copyto!
     ╎    ╎    ╎    ╎    ╎    ╎   9      …actarray.jl:1240; pointer
    9╎    ╎    ╎    ╎    ╎    ╎    9      …/pointer.jl:282; +
     ╎    ╎    ╎    ╎    ╎    ╎  64     …se/array.jl:337; unsafe_copyto!
   46╎    ╎    ╎    ╎    ╎    ╎   64     …ase/cmem.jl:26; memmove
     ╎    ╎    ╎    ╎    ╎    ╎    18     …sentials.jl:543; cconvert
     ╎    ╎    ╎    ╎    ╎    ╎     18     …e/number.jl:7; convert
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 18     …ase/boot.jl:789; UInt64
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  18     …ase/boot.jl:759; toUInt64
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   18     …ase/boot.jl:648; check_top_bit
   18╎    ╎    ╎    ╎    ╎    ╎    ╎    18     …ase/boot.jl:638; is_top_bit_set
    2╎    ╎    ╎    ╎    ╎    2      …sentials.jl:13; getindex
   18╎    ╎    ╎    ╎    ╎    626    …/generic.jl:596; norm
   13╎    ╎    ╎    ╎    ╎     13     …/generic.jl:0; generic_norm2(x::Vector{F…
   27╎    ╎    ╎    ╎    ╎     27     …/generic.jl:595; norm(itr::Vector{Float6…
     ╎    ╎    ╎    ╎    ╎     5      …/generic.jl:596; norm(itr::Vector{Float6…
     ╎    ╎    ╎    ╎    ╎    ╎ 5      …actarray.jl:1220; isempty
    5╎    ╎    ╎    ╎    ╎    ╎  5      …sentials.jl:10; length
   10╎    ╎    ╎    ╎    ╎     10     …/generic.jl:597; norm(itr::Vector{Float6…
     ╎    ╎    ╎    ╎    ╎     553    …/generic.jl:598; norm(itr::Vector{Float6…
    8╎    ╎    ╎    ╎    ╎    ╎ 553    …rc/dense.jl:106; norm2
   25╎    ╎    ╎    ╎    ╎    ╎  25     …se/array.jl:0; generic_norm2(x::Vector…
    8╎    ╎    ╎    ╎    ╎    ╎  8      …e/reduce.jl:0; _mapreduce(f::typeof(Li…
    6╎    ╎    ╎    ╎    ╎    ╎  6      …e/reduce.jl:428; _mapreduce(f::typeof(…
   21╎    ╎    ╎    ╎    ╎    ╎  21     …/generic.jl:0; generic_norm2(x::Vector…
   19╎    ╎    ╎    ╎    ╎    ╎  19     …/generic.jl:462; generic_norm2(x::Vect…
     ╎    ╎    ╎    ╎    ╎    ╎  300    …/generic.jl:463; generic_norm2(x::Vect…
     ╎    ╎    ╎    ╎    ╎    ╎   300    …/generic.jl:527; normInf
     ╎    ╎    ╎    ╎    ╎    ╎    300    …/generic.jl:453; generic_normInf
     ╎    ╎    ╎    ╎    ╎    ╎     300    …educedim.jl:357; mapreduce
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 300    …educedim.jl:357; #mapreduce#821
   13╎    ╎    ╎    ╎    ╎    ╎    ╎  300    …educedim.jl:365; _mapreduce_dim
   20╎    ╎    ╎    ╎    ╎    ╎    ╎   20     …e/reduce.jl:0; _mapreduce(f::typ…
   17╎    ╎    ╎    ╎    ╎    ╎    ╎   17     …e/reduce.jl:428; _mapreduce(f::t…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   6      …e/reduce.jl:429; _mapreduce(f::t…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    6      …/indices.jl:486; LinearIndices
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     6      …actarray.jl:98; axes
    6╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 6      …se/array.jl:191; size
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   2      …e/reduce.jl:436; _mapreduce(f::t…
    2╎    ╎    ╎    ╎    ╎    ╎    ╎    2      @Base/int.jl:83; <
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   133    …e/reduce.jl:440; _mapreduce(f::t…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    37     …ase/math.jl:908; max
   37╎    ╎    ╎    ╎    ╎    ╎    ╎     37     …se/float.jl:410; -
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    65     …ase/math.jl:909; max
   25╎    ╎    ╎    ╎    ╎    ╎    ╎     25     …sentials.jl:647; ifelse
   40╎    ╎    ╎    ╎    ╎    ╎    ╎     40     …oatfuncs.jl:15; signbit
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    13     …ase/math.jl:911; max
   13╎    ╎    ╎    ╎    ╎    ╎    ╎     13     …sentials.jl:647; ifelse
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    18     …/generic.jl:639; norm
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     18     …/generic.jl:639; norm
   18╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 18     …se/float.jl:610; abs
   37╎    ╎    ╎    ╎    ╎    ╎    ╎   41     …e/reduce.jl:441; _mapreduce(f::t…
    4╎    ╎    ╎    ╎    ╎    ╎    ╎    4      @Base/int.jl:83; <
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   3      …e/reduce.jl:442; _mapreduce(f::t…
    3╎    ╎    ╎    ╎    ╎    ╎    ╎    3      …sentials.jl:13; getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   65     …e/reduce.jl:443; _mapreduce(f::t…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    9      …ase/math.jl:908; max
    9╎    ╎    ╎    ╎    ╎    ╎    ╎     9      …se/float.jl:410; -
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    42     …ase/math.jl:909; max
   28╎    ╎    ╎    ╎    ╎    ╎    ╎     28     …sentials.jl:647; ifelse
   14╎    ╎    ╎    ╎    ╎    ╎    ╎     14     …oatfuncs.jl:15; signbit
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    14     …ase/math.jl:911; max
   14╎    ╎    ╎    ╎    ╎    ╎    ╎     14     …sentials.jl:647; ifelse
    1╎    ╎    ╎    ╎    ╎    ╎  29     …/generic.jl:464; generic_norm2(x::Vect…
   12╎    ╎    ╎    ╎    ╎    ╎   28     …se/float.jl:635; isinf
   16╎    ╎    ╎    ╎    ╎    ╎    16     …se/float.jl:610; abs
     ╎    ╎    ╎    ╎    ╎    ╎  8      …/generic.jl:465; generic_norm2(x::Vect…
     ╎    ╎    ╎    ╎    ╎    ╎   8      …se/array.jl:945; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    8      …se/array.jl:945; iterate
    8╎    ╎    ╎    ╎    ╎    ╎     8      …sentials.jl:13; getindex
   13╎    ╎    ╎    ╎    ╎    ╎  73     …/generic.jl:467; generic_norm2(x::Vect…
     ╎    ╎    ╎    ╎    ╎    ╎   48     …se/float.jl:623; isfinite
   42╎    ╎    ╎    ╎    ╎    ╎    42     …se/float.jl:410; -
     ╎    ╎    ╎    ╎    ╎    ╎    6      …se/float.jl:620; isnan
    6╎    ╎    ╎    ╎    ╎    ╎     6      …se/float.jl:535; !=
     ╎    ╎    ╎    ╎    ╎    ╎   12     …perators.jl:587; *
   12╎    ╎    ╎    ╎    ╎    ╎    12     …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎  13     …/generic.jl:468; generic_norm2(x::Vect…
     ╎    ╎    ╎    ╎    ╎    ╎   13     …/generic.jl:459; norm_sqr
     ╎    ╎    ╎    ╎    ╎    ╎    13     …e/number.jl:189; abs2
   13╎    ╎    ╎    ╎    ╎    ╎     13     …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎  30     …/generic.jl:470; generic_norm2(x::Vect…
    9╎    ╎    ╎    ╎    ╎    ╎   30     …se/array.jl:945; iterate
    7╎    ╎    ╎    ╎    ╎    ╎    7      …sentials.jl:13; getindex
   14╎    ╎    ╎    ╎    ╎    ╎    14     @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎    ╎  13     …/generic.jl:473; generic_norm2(x::Vect…
   13╎    ╎    ╎    ╎    ╎    ╎   13     …se/float.jl:409; +
     ╎    ╎    ╎    ╎    ╎    68     …c/matmul.jl:15; dot
     ╎    ╎    ╎    ╎    ╎     8      …src/blas.jl:393; dot
    8╎    ╎    ╎    ╎    ╎    ╎ 8      …sentials.jl:10; length
     ╎    ╎    ╎    ╎    ╎     1      …src/blas.jl:394; dot
    1╎    ╎    ╎    ╎    ╎    ╎ 1      …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎     59     …src/blas.jl:395; dot
   55╎    ╎    ╎    ╎    ╎    ╎ 55     …src/blas.jl:345; dot
     ╎    ╎    ╎    ╎    ╎    ╎ 4      …src/blas.jl:177; vec_pointer_stride
     ╎    ╎    ╎    ╎    ╎    ╎  4      …src/blas.jl:177; vec_pointer_stride
     ╎    ╎    ╎    ╎    ╎    ╎   4      …actarray.jl:1237; pointer
    4╎    ╎    ╎    ╎    ╎    ╎    4      …/pointer.jl:65; unsafe_convert
     ╎    ╎    ╎    ╎    ╎   198    …_physics.jl:72; projected!(A::Matrix{Vecto…
    9╎    ╎    ╎    ╎    ╎    9      …se/float.jl:411; *
    3╎    ╎    ╎    ╎    ╎    189    …/generic.jl:596; norm
    3╎    ╎    ╎    ╎    ╎     3      …/generic.jl:0; generic_norm2(x::Vector{F…
    5╎    ╎    ╎    ╎    ╎     5      …/generic.jl:595; norm(itr::Vector{Float6…
     ╎    ╎    ╎    ╎    ╎     5      …/generic.jl:596; norm(itr::Vector{Float6…
     ╎    ╎    ╎    ╎    ╎    ╎ 5      …actarray.jl:1220; isempty
    5╎    ╎    ╎    ╎    ╎    ╎  5      …sentials.jl:10; length
    7╎    ╎    ╎    ╎    ╎     7      …/generic.jl:597; norm(itr::Vector{Float6…
     ╎    ╎    ╎    ╎    ╎     166    …/generic.jl:598; norm(itr::Vector{Float6…
     ╎    ╎    ╎    ╎    ╎    ╎ 166    …rc/dense.jl:106; norm2
   14╎    ╎    ╎    ╎    ╎    ╎  14     …se/array.jl:0; generic_norm2(x::Vector…
    5╎    ╎    ╎    ╎    ╎    ╎  5      …e/reduce.jl:428; _mapreduce(f::typeof(…
   13╎    ╎    ╎    ╎    ╎    ╎  13     …/generic.jl:0; generic_norm2(x::Vector…
   18╎    ╎    ╎    ╎    ╎    ╎  18     …/generic.jl:462; generic_norm2(x::Vect…
     ╎    ╎    ╎    ╎    ╎    ╎  64     …/generic.jl:463; generic_norm2(x::Vect…
     ╎    ╎    ╎    ╎    ╎    ╎   64     …/generic.jl:527; normInf
     ╎    ╎    ╎    ╎    ╎    ╎    64     …/generic.jl:453; generic_normInf
     ╎    ╎    ╎    ╎    ╎    ╎     64     …educedim.jl:357; mapreduce
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 64     …educedim.jl:357; #mapreduce#821
    6╎    ╎    ╎    ╎    ╎    ╎    ╎  64     …educedim.jl:365; _mapreduce_dim
    3╎    ╎    ╎    ╎    ╎    ╎    ╎   3      …e/reduce.jl:0; _mapreduce(f::typ…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   2      …e/reduce.jl:429; _mapreduce(f::t…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    2      …/indices.jl:486; LinearIndices
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     2      …actarray.jl:98; axes
    2╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 2      …se/array.jl:191; size
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   23     …e/reduce.jl:440; _mapreduce(f::t…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    6      …ase/math.jl:908; max
    6╎    ╎    ╎    ╎    ╎    ╎    ╎     6      …se/float.jl:410; -
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    8      …ase/math.jl:909; max
    8╎    ╎    ╎    ╎    ╎    ╎    ╎     8      …sentials.jl:647; ifelse
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    9      …ase/math.jl:911; max
    9╎    ╎    ╎    ╎    ╎    ╎    ╎     9      …sentials.jl:647; ifelse
   13╎    ╎    ╎    ╎    ╎    ╎    ╎   13     …e/reduce.jl:441; _mapreduce(f::t…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   17     …e/reduce.jl:443; _mapreduce(f::t…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    5      …ase/math.jl:909; max
    5╎    ╎    ╎    ╎    ╎    ╎    ╎     5      …sentials.jl:647; ifelse
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    7      …ase/math.jl:911; max
    7╎    ╎    ╎    ╎    ╎    ╎    ╎     7      …sentials.jl:647; ifelse
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    5      …/generic.jl:639; norm
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     5      …/generic.jl:639; norm
    5╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 5      …se/float.jl:610; abs
     ╎    ╎    ╎    ╎    ╎    ╎  14     …/generic.jl:464; generic_norm2(x::Vect…
    8╎    ╎    ╎    ╎    ╎    ╎   14     …se/float.jl:635; isinf
    6╎    ╎    ╎    ╎    ╎    ╎    6      …se/float.jl:610; abs
     ╎    ╎    ╎    ╎    ╎    ╎  5      …/generic.jl:465; generic_norm2(x::Vect…
     ╎    ╎    ╎    ╎    ╎    ╎   5      …se/array.jl:945; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    5      …se/array.jl:945; iterate
    5╎    ╎    ╎    ╎    ╎    ╎     5      …sentials.jl:13; getindex
    1╎    ╎    ╎    ╎    ╎    ╎  8      …/generic.jl:467; generic_norm2(x::Vect…
     ╎    ╎    ╎    ╎    ╎    ╎   2      …se/float.jl:623; isfinite
    1╎    ╎    ╎    ╎    ╎    ╎    1      …se/float.jl:410; -
     ╎    ╎    ╎    ╎    ╎    ╎    1      …se/float.jl:620; isnan
    1╎    ╎    ╎    ╎    ╎    ╎     1      …se/float.jl:535; !=
     ╎    ╎    ╎    ╎    ╎    ╎   5      …perators.jl:587; *
    5╎    ╎    ╎    ╎    ╎    ╎    5      …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎  3      …/generic.jl:468; generic_norm2(x::Vect…
     ╎    ╎    ╎    ╎    ╎    ╎   3      …/generic.jl:459; norm_sqr
     ╎    ╎    ╎    ╎    ╎    ╎    3      …e/number.jl:189; abs2
    3╎    ╎    ╎    ╎    ╎    ╎     3      …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎  15     …/generic.jl:470; generic_norm2(x::Vect…
    5╎    ╎    ╎    ╎    ╎    ╎   15     …se/array.jl:945; iterate
    5╎    ╎    ╎    ╎    ╎    ╎    5      …sentials.jl:13; getindex
    5╎    ╎    ╎    ╎    ╎    ╎    5      @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎    ╎  7      …/generic.jl:473; generic_norm2(x::Vect…
    7╎    ╎    ╎    ╎    ╎    ╎   7      …se/float.jl:409; +
    4╎    ╎    ╎    ╎    ╎   4      …_physics.jl:73; projected!(A::Matrix{Vecto…
    1╎    ╎    ╎    ╎    ╎   1      …_physics.jl:74; projected!(A::Matrix{Vecto…
     ╎    ╎    ╎    ╎    ╎  4      …pse_comp.jl:91; eclipse_compute_quantities!…
     ╎    ╎    ╎    ╎    ╎   4      …roadcast.jl:911; materialize!
     ╎    ╎    ╎    ╎    ╎    4      …roadcast.jl:914; materialize!
     ╎    ╎    ╎    ╎    ╎     4      …roadcast.jl:956; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎ 4      …roadcast.jl:1003; copyto!
    1╎    ╎    ╎    ╎    ╎    ╎  1      …simdloop.jl:0; macro expansion
    1╎    ╎    ╎    ╎    ╎    ╎  1      …simdloop.jl:75; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎  2      …simdloop.jl:77; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎   2      …roadcast.jl:1004; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎    1      …roadcast.jl:636; getindex
     ╎    ╎    ╎    ╎    ╎    ╎     1      …roadcast.jl:682; _broadcast_getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …roadcast.jl:709; _broadcast_getind…
    1╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎    1      …ensional.jl:698; setindex!
    1╎    ╎    ╎    ╎    ╎    ╎     1      …se/array.jl:1024; setindex!
     ╎    ╎    ╎    ╎    ╎  40     …pse_comp.jl:92; eclipse_compute_quantities!…
     ╎    ╎    ╎    ╎    ╎   40     …roadcast.jl:903; materialize
     ╎    ╎    ╎    ╎    ╎    40     …roadcast.jl:928; copy
     ╎    ╎    ╎    ╎    ╎     20     …roadcast.jl:956; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎ 19     …roadcast.jl:1003; copyto!
    1╎    ╎    ╎    ╎    ╎    ╎  1      @Base/int.jl:0; macro expansion
    1╎    ╎    ╎    ╎    ╎    ╎  1      …simdloop.jl:0; macro expansion
    1╎    ╎    ╎    ╎    ╎    ╎  1      …simdloop.jl:75; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎  16     …simdloop.jl:77; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎   16     …roadcast.jl:1004; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎    7      …roadcast.jl:636; getindex
     ╎    ╎    ╎    ╎    ╎    ╎     1      …roadcast.jl:681; _broadcast_getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …roadcast.jl:705; _getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …roadcast.jl:675; _broadcast_getin…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …ensional.jl:696; getindex
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …sentials.jl:14; getindex
     ╎    ╎    ╎    ╎    ╎    ╎     6      …roadcast.jl:682; _broadcast_getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 6      …roadcast.jl:709; _broadcast_getind…
    6╎    ╎    ╎    ╎    ╎    ╎    ╎  6      …se/float.jl:412; /
     ╎    ╎    ╎    ╎    ╎    ╎    9      …ensional.jl:698; setindex!
    9╎    ╎    ╎    ╎    ╎    ╎     9      …se/array.jl:1024; setindex!
    1╎    ╎    ╎    ╎    ╎    ╎ 1      …simdloop.jl:0; copyto!
     ╎    ╎    ╎    ╎    ╎     20     …roadcast.jl:223; similar
     ╎    ╎    ╎    ╎    ╎    ╎ 20     …roadcast.jl:224; similar
     ╎    ╎    ╎    ╎    ╎    ╎  20     …actarray.jl:876; similar
     ╎    ╎    ╎    ╎    ╎    ╎   20     …actarray.jl:877; similar
     ╎    ╎    ╎    ╎    ╎    ╎    20     …ase/boot.jl:494; Array
     ╎    ╎    ╎    ╎    ╎    ╎     20     …ase/boot.jl:487; Array
   20╎    ╎    ╎    ╎    ╎    ╎    ╎ 20     …ase/boot.jl:479; Array
     ╎    ╎    ╎    ╎    ╎  40     …pse_comp.jl:95; eclipse_compute_quantities!…
     ╎    ╎    ╎    ╎    ╎   1      …se/array.jl:973; getindex
    1╎    ╎    ╎    ╎    ╎    1      …actarray.jl:700; checkbounds
     ╎    ╎    ╎    ╎    ╎   39     …se/array.jl:975; getindex
     ╎    ╎    ╎    ╎    ╎    39     …actarray.jl:831; similar
     ╎    ╎    ╎    ╎    ╎     39     …se/array.jl:420; similar
     ╎    ╎    ╎    ╎    ╎    ╎ 39     …ase/boot.jl:486; Array
   39╎    ╎    ╎    ╎    ╎    ╎  39     …ase/boot.jl:477; Array
     ╎    ╎    ╎    ╎    ╎  12     …pse_comp.jl:97; eclipse_compute_quantities!…
    5╎    ╎    ╎    ╎    ╎   5      …sentials.jl:0; getindex
     ╎    ╎    ╎    ╎    ╎   7      …subarray.jl:184; view
    7╎    ╎    ╎    ╎    ╎    7      …actarray.jl:702; checkbounds
     ╎    ╎    ╎    ╎    ╎  673    …pse_comp.jl:98; eclipse_compute_quantities!…
    7╎    ╎    ╎    ╎    ╎   7      …se/float.jl:411; *
   22╎    ╎    ╎    ╎    ╎   582    …/generic.jl:596; norm
    4╎    ╎    ╎    ╎    ╎    4      …/generic.jl:0; generic_norm2(x::Vector{Fl…
   10╎    ╎    ╎    ╎    ╎    10     …/generic.jl:595; norm(itr::Vector{Float64…
    6╎    ╎    ╎    ╎    ╎    6      …/generic.jl:597; norm(itr::Vector{Float64…
     ╎    ╎    ╎    ╎    ╎    328    …/generic.jl:598; norm
   16╎    ╎    ╎    ╎    ╎     328    …rc/dense.jl:106; norm2
    6╎    ╎    ╎    ╎    ╎    ╎ 6      …e/reduce.jl:0; _mapreduce(f::typeof(Lin…
   10╎    ╎    ╎    ╎    ╎    ╎ 10     …e/reduce.jl:428; _mapreduce(f::typeof(L…
   26╎    ╎    ╎    ╎    ╎    ╎ 26     …/generic.jl:0; generic_norm2(x::SubArra…
   13╎    ╎    ╎    ╎    ╎    ╎ 13     …/generic.jl:462; generic_norm2(x::SubAr…
     ╎    ╎    ╎    ╎    ╎    ╎ 134    …/generic.jl:463; generic_norm2(x::SubAr…
     ╎    ╎    ╎    ╎    ╎    ╎  134    …/generic.jl:527; normInf
     ╎    ╎    ╎    ╎    ╎    ╎   134    …/generic.jl:453; generic_normInf
     ╎    ╎    ╎    ╎    ╎    ╎    134    …educedim.jl:357; mapreduce
     ╎    ╎    ╎    ╎    ╎    ╎     134    …educedim.jl:357; #mapreduce#821
    4╎    ╎    ╎    ╎    ╎    ╎    ╎ 134    …educedim.jl:365; _mapreduce_dim
   25╎    ╎    ╎    ╎    ╎    ╎    ╎  25     …e/reduce.jl:0; _mapreduce(f::type…
   19╎    ╎    ╎    ╎    ╎    ╎    ╎  19     …e/reduce.jl:428; _mapreduce(f::ty…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  5      …e/reduce.jl:429; _mapreduce(f::ty…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   5      …/indices.jl:486; LinearIndices
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    5      …subarray.jl:490; axes
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     5      …subarray.jl:495; _indices_sub
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 5      …se/range.jl:706; axes
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  5      …se/range.jl:761; length
    5╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   5      @Base/int.jl:86; -
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …e/reduce.jl:431; _mapreduce(f::ty…
    1╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  5      …e/reduce.jl:433; _mapreduce(f::ty…
    5╎    ╎    ╎    ╎    ╎    ╎    ╎   5      …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  8      …e/reduce.jl:438; _mapreduce(f::ty…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   8      …subarray.jl:323; getindex
    8╎    ╎    ╎    ╎    ╎    ╎    ╎    8      …ase/Base.jl:37; getproperty
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  20     …e/reduce.jl:440; _mapreduce(f::ty…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …ase/math.jl:908; max
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …se/float.jl:410; -
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   7      …ase/math.jl:909; max
    7╎    ╎    ╎    ╎    ╎    ╎    ╎    7      …sentials.jl:647; ifelse
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   10     …ase/math.jl:911; max
   10╎    ╎    ╎    ╎    ╎    ╎    ╎    10     …sentials.jl:647; ifelse
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   2      …/generic.jl:639; norm
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    2      …/generic.jl:639; norm
    2╎    ╎    ╎    ╎    ╎    ╎    ╎     2      …se/float.jl:610; abs
   33╎    ╎    ╎    ╎    ╎    ╎    ╎  33     …e/reduce.jl:441; _mapreduce(f::ty…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …e/reduce.jl:442; _mapreduce(f::ty…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …subarray.jl:323; getindex
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …sentials.jl:13; getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  13     …e/reduce.jl:443; _mapreduce(f::ty…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   7      …ase/math.jl:909; max
    7╎    ╎    ╎    ╎    ╎    ╎    ╎    7      …sentials.jl:647; ifelse
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   4      …ase/math.jl:911; max
    4╎    ╎    ╎    ╎    ╎    ╎    ╎    4      …sentials.jl:647; ifelse
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   2      …/generic.jl:639; norm
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    2      …/generic.jl:639; norm
    2╎    ╎    ╎    ╎    ╎    ╎    ╎     2      …se/float.jl:610; abs
     ╎    ╎    ╎    ╎    ╎    ╎ 13     …/generic.jl:464; generic_norm2(x::SubAr…
    8╎    ╎    ╎    ╎    ╎    ╎  12     …se/float.jl:635; isinf
    4╎    ╎    ╎    ╎    ╎    ╎   4      …se/float.jl:610; abs
     ╎    ╎    ╎    ╎    ╎    ╎  1      …e/number.jl:42; iszero
    1╎    ╎    ╎    ╎    ╎    ╎   1      …se/float.jl:534; ==
     ╎    ╎    ╎    ╎    ╎    ╎ 28     …/generic.jl:465; generic_norm2(x::SubAr…
    1╎    ╎    ╎    ╎    ╎    ╎  1      …actarray.jl:0; iterate
     ╎    ╎    ╎    ╎    ╎    ╎  27     …actarray.jl:1215; iterate
     ╎    ╎    ╎    ╎    ╎    ╎   2      …actarray.jl:321; eachindex
     ╎    ╎    ╎    ╎    ╎    ╎    2      …actarray.jl:137; axes1
     ╎    ╎    ╎    ╎    ╎    ╎     2      …subarray.jl:490; axes
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 2      …subarray.jl:495; _indices_sub
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  2      …se/range.jl:706; axes
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   2      …se/range.jl:761; length
    2╎    ╎    ╎    ╎    ╎    ╎    ╎    2      @Base/int.jl:86; -
    7╎    ╎    ╎    ╎    ╎    ╎   7      …actarray.jl:0; ifelse
     ╎    ╎    ╎    ╎    ╎    ╎   9      …actarray.jl:1215; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    9      …se/range.jl:897; iterate
     ╎    ╎    ╎    ╎    ╎    ╎     9      …se/range.jl:672; isempty
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 9      …perators.jl:378; >
    9╎    ╎    ╎    ╎    ╎    ╎    ╎  9      @Base/int.jl:83; <
     ╎    ╎    ╎    ╎    ╎    ╎   9      …actarray.jl:1217; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    9      …subarray.jl:323; getindex
    9╎    ╎    ╎    ╎    ╎    ╎     9      …sentials.jl:13; getindex
    9╎    ╎    ╎    ╎    ╎    ╎ 35     …/generic.jl:467; generic_norm2(x::SubAr…
     ╎    ╎    ╎    ╎    ╎    ╎  17     …se/float.jl:623; isfinite
    2╎    ╎    ╎    ╎    ╎    ╎   2      …se/float.jl:410; -
     ╎    ╎    ╎    ╎    ╎    ╎   15     …se/float.jl:620; isnan
   15╎    ╎    ╎    ╎    ╎    ╎    15     …se/float.jl:535; !=
     ╎    ╎    ╎    ╎    ╎    ╎  9      …perators.jl:587; *
    7╎    ╎    ╎    ╎    ╎    ╎   7      …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎   2      …romotion.jl:423; *
    2╎    ╎    ╎    ╎    ╎    ╎    2      …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎ 11     …/generic.jl:468; generic_norm2(x::SubAr…
     ╎    ╎    ╎    ╎    ╎    ╎  11     …/generic.jl:459; norm_sqr
     ╎    ╎    ╎    ╎    ╎    ╎   11     …e/number.jl:189; abs2
   11╎    ╎    ╎    ╎    ╎    ╎    11     …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎ 27     …/generic.jl:470; generic_norm2(x::SubAr…
     ╎    ╎    ╎    ╎    ╎    ╎  12     …actarray.jl:1215; iterate
     ╎    ╎    ╎    ╎    ╎    ╎   12     …se/range.jl:901; iterate
   12╎    ╎    ╎    ╎    ╎    ╎    12     …romotion.jl:521; ==
    6╎    ╎    ╎    ╎    ╎    ╎  6      …actarray.jl:1216; iterate
     ╎    ╎    ╎    ╎    ╎    ╎  9      …actarray.jl:1217; iterate
     ╎    ╎    ╎    ╎    ╎    ╎   9      …subarray.jl:323; getindex
    9╎    ╎    ╎    ╎    ╎    ╎    9      …sentials.jl:13; getindex
     ╎    ╎    ╎    ╎    ╎    ╎ 5      …/generic.jl:473; generic_norm2(x::SubAr…
    5╎    ╎    ╎    ╎    ╎    ╎  5      …se/float.jl:409; +
     ╎    ╎    ╎    ╎    ╎    ╎ 4      …/generic.jl:476; generic_norm2(x::SubAr…
    4╎    ╎    ╎    ╎    ╎    ╎  4      …ase/math.jl:686; sqrt
     ╎    ╎    ╎    ╎    ╎    212    …/generic.jl:598; norm(itr::Vector{Float64…
     ╎    ╎    ╎    ╎    ╎     212    …rc/dense.jl:106; norm2
   10╎    ╎    ╎    ╎    ╎    ╎ 10     …se/array.jl:0; generic_norm2(x::Vector{…
    6╎    ╎    ╎    ╎    ╎    ╎ 6      …e/reduce.jl:428; _mapreduce(f::typeof(L…
   16╎    ╎    ╎    ╎    ╎    ╎ 16     …/generic.jl:0; generic_norm2(x::Vector{…
   16╎    ╎    ╎    ╎    ╎    ╎ 16     …/generic.jl:462; generic_norm2(x::Vecto…
     ╎    ╎    ╎    ╎    ╎    ╎ 95     …/generic.jl:463; generic_norm2(x::Vecto…
     ╎    ╎    ╎    ╎    ╎    ╎  95     …/generic.jl:527; normInf
     ╎    ╎    ╎    ╎    ╎    ╎   95     …/generic.jl:453; generic_normInf
     ╎    ╎    ╎    ╎    ╎    ╎    95     …educedim.jl:357; mapreduce
     ╎    ╎    ╎    ╎    ╎    ╎     95     …educedim.jl:357; #mapreduce#821
   11╎    ╎    ╎    ╎    ╎    ╎    ╎ 95     …educedim.jl:365; _mapreduce_dim
   12╎    ╎    ╎    ╎    ╎    ╎    ╎  12     …e/reduce.jl:0; _mapreduce(f::type…
    7╎    ╎    ╎    ╎    ╎    ╎    ╎  7      …e/reduce.jl:428; _mapreduce(f::ty…
    7╎    ╎    ╎    ╎    ╎    ╎    ╎  7      …e/reduce.jl:431; _mapreduce(f::ty…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  3      …e/reduce.jl:438; _mapreduce(f::ty…
    3╎    ╎    ╎    ╎    ╎    ╎    ╎   3      …sentials.jl:13; getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  18     …e/reduce.jl:440; _mapreduce(f::ty…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   5      …ase/math.jl:908; max
    5╎    ╎    ╎    ╎    ╎    ╎    ╎    5      …se/float.jl:410; -
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   12     …ase/math.jl:909; max
    6╎    ╎    ╎    ╎    ╎    ╎    ╎    6      …sentials.jl:647; ifelse
    6╎    ╎    ╎    ╎    ╎    ╎    ╎    6      …oatfuncs.jl:15; signbit
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …ase/math.jl:911; max
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …sentials.jl:647; ifelse
   17╎    ╎    ╎    ╎    ╎    ╎    ╎  20     …e/reduce.jl:441; _mapreduce(f::ty…
    3╎    ╎    ╎    ╎    ╎    ╎    ╎   3      @Base/int.jl:83; <
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  17     …e/reduce.jl:443; _mapreduce(f::ty…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   6      …ase/math.jl:908; max
    6╎    ╎    ╎    ╎    ╎    ╎    ╎    6      …se/float.jl:410; -
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   9      …ase/math.jl:909; max
    9╎    ╎    ╎    ╎    ╎    ╎    ╎    9      …sentials.jl:647; ifelse
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   2      …ase/math.jl:911; max
    2╎    ╎    ╎    ╎    ╎    ╎    ╎    2      …sentials.jl:647; ifelse
    1╎    ╎    ╎    ╎    ╎    ╎ 18     …/generic.jl:464; generic_norm2(x::Vecto…
   10╎    ╎    ╎    ╎    ╎    ╎  17     …se/float.jl:635; isinf
    7╎    ╎    ╎    ╎    ╎    ╎   7      …se/float.jl:610; abs
     ╎    ╎    ╎    ╎    ╎    ╎ 9      …/generic.jl:465; generic_norm2(x::Vecto…
     ╎    ╎    ╎    ╎    ╎    ╎  9      …se/array.jl:945; iterate
     ╎    ╎    ╎    ╎    ╎    ╎   9      …se/array.jl:945; iterate
    8╎    ╎    ╎    ╎    ╎    ╎    8      …sentials.jl:13; getindex
    1╎    ╎    ╎    ╎    ╎    ╎    1      …sentials.jl:10; length
    3╎    ╎    ╎    ╎    ╎    ╎ 13     …/generic.jl:467; generic_norm2(x::Vecto…
     ╎    ╎    ╎    ╎    ╎    ╎  3      …se/float.jl:623; isfinite
    3╎    ╎    ╎    ╎    ╎    ╎   3      …se/float.jl:410; -
     ╎    ╎    ╎    ╎    ╎    ╎  7      …perators.jl:587; *
    6╎    ╎    ╎    ╎    ╎    ╎   6      …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎   1      …romotion.jl:423; *
     ╎    ╎    ╎    ╎    ╎    ╎    1      …romotion.jl:393; promote
     ╎    ╎    ╎    ╎    ╎    ╎     1      …romotion.jl:370; _promote
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …e/number.jl:7; convert
    1╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …se/float.jl:159; Float64
     ╎    ╎    ╎    ╎    ╎    ╎ 3      …/generic.jl:468; generic_norm2(x::Vecto…
     ╎    ╎    ╎    ╎    ╎    ╎  3      …/generic.jl:459; norm_sqr
     ╎    ╎    ╎    ╎    ╎    ╎   3      …e/number.jl:189; abs2
    3╎    ╎    ╎    ╎    ╎    ╎    3      …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎ 22     …/generic.jl:470; generic_norm2(x::Vecto…
   13╎    ╎    ╎    ╎    ╎    ╎  22     …se/array.jl:945; iterate
    5╎    ╎    ╎    ╎    ╎    ╎   5      …sentials.jl:13; getindex
    4╎    ╎    ╎    ╎    ╎    ╎   4      @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎    ╎ 4      …/generic.jl:473; generic_norm2(x::Vecto…
    4╎    ╎    ╎    ╎    ╎    ╎  4      …se/float.jl:409; +
     ╎    ╎    ╎    ╎    ╎   84     …c/matmul.jl:15; dot
     ╎    ╎    ╎    ╎    ╎    1      …src/blas.jl:393; dot
    1╎    ╎    ╎    ╎    ╎     1      …sentials.jl:10; length
     ╎    ╎    ╎    ╎    ╎    83     …src/blas.jl:395; dot
   74╎    ╎    ╎    ╎    ╎     74     …src/blas.jl:345; dot
     ╎    ╎    ╎    ╎    ╎     9      …src/blas.jl:177; vec_pointer_stride
     ╎    ╎    ╎    ╎    ╎    ╎ 9      …src/blas.jl:177; vec_pointer_stride
     ╎    ╎    ╎    ╎    ╎    ╎  9      …actarray.jl:1237; pointer
    8╎    ╎    ╎    ╎    ╎    ╎   8      …/pointer.jl:65; unsafe_convert
     ╎    ╎    ╎    ╎    ╎    ╎   1      …subarray.jl:472; unsafe_convert
    1╎    ╎    ╎    ╎    ╎    ╎    1      …/pointer.jl:65; unsafe_convert
     ╎    ╎    ╎    ╎    ╎  257    …pse_comp.jl:99; eclipse_compute_quantities!…
    1╎    ╎    ╎    ╎    ╎   1      …se/array.jl:1021; setindex!
   10╎    ╎    ╎    ╎    ╎   10     …se/float.jl:411; *
    4╎    ╎    ╎    ╎    ╎   246    …/generic.jl:596; norm
    7╎    ╎    ╎    ╎    ╎    7      …/generic.jl:0; generic_norm2(x::Vector{Fl…
    4╎    ╎    ╎    ╎    ╎    4      …/generic.jl:595; norm(itr::Vector{Float64…
     ╎    ╎    ╎    ╎    ╎    6      …/generic.jl:596; norm(itr::Vector{Float64…
     ╎    ╎    ╎    ╎    ╎     6      …actarray.jl:1220; isempty
    4╎    ╎    ╎    ╎    ╎    ╎ 4      …sentials.jl:10; length
    2╎    ╎    ╎    ╎    ╎    ╎ 2      …romotion.jl:521; ==
    6╎    ╎    ╎    ╎    ╎    6      …/generic.jl:597; norm(itr::Vector{Float64…
     ╎    ╎    ╎    ╎    ╎    219    …/generic.jl:598; norm(itr::Vector{Float64…
     ╎    ╎    ╎    ╎    ╎     219    …rc/dense.jl:106; norm2
   17╎    ╎    ╎    ╎    ╎    ╎ 17     …se/array.jl:0; generic_norm2(x::Vector{…
    1╎    ╎    ╎    ╎    ╎    ╎ 1      @Base/int.jl:83; <
    8╎    ╎    ╎    ╎    ╎    ╎ 8      …e/reduce.jl:428; _mapreduce(f::typeof(L…
   35╎    ╎    ╎    ╎    ╎    ╎ 35     …/generic.jl:0; generic_norm2(x::Vector{…
   13╎    ╎    ╎    ╎    ╎    ╎ 13     …/generic.jl:462; generic_norm2(x::Vecto…
     ╎    ╎    ╎    ╎    ╎    ╎ 74     …/generic.jl:463; generic_norm2(x::Vecto…
     ╎    ╎    ╎    ╎    ╎    ╎  74     …/generic.jl:527; normInf
     ╎    ╎    ╎    ╎    ╎    ╎   74     …/generic.jl:453; generic_normInf
     ╎    ╎    ╎    ╎    ╎    ╎    74     …educedim.jl:357; mapreduce
     ╎    ╎    ╎    ╎    ╎    ╎     74     …educedim.jl:357; #mapreduce#821
    5╎    ╎    ╎    ╎    ╎    ╎    ╎ 74     …educedim.jl:365; _mapreduce_dim
    8╎    ╎    ╎    ╎    ╎    ╎    ╎  8      …e/reduce.jl:0; _mapreduce(f::type…
    5╎    ╎    ╎    ╎    ╎    ╎    ╎  5      …e/reduce.jl:428; _mapreduce(f::ty…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  4      …e/reduce.jl:429; _mapreduce(f::ty…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   4      …/indices.jl:486; LinearIndices
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    4      …actarray.jl:98; axes
    4╎    ╎    ╎    ╎    ╎    ╎    ╎     4      …se/array.jl:191; size
    1╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …e/reduce.jl:431; _mapreduce(f::ty…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  18     …e/reduce.jl:440; _mapreduce(f::ty…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   7      …ase/math.jl:908; max
    7╎    ╎    ╎    ╎    ╎    ╎    ╎    7      …se/float.jl:410; -
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   3      …ase/math.jl:909; max
    3╎    ╎    ╎    ╎    ╎    ╎    ╎    3      …sentials.jl:647; ifelse
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   8      …ase/math.jl:911; max
    8╎    ╎    ╎    ╎    ╎    ╎    ╎    8      …sentials.jl:647; ifelse
   16╎    ╎    ╎    ╎    ╎    ╎    ╎  16     …e/reduce.jl:441; _mapreduce(f::ty…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  17     …e/reduce.jl:443; _mapreduce(f::ty…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   7      …ase/math.jl:909; max
    7╎    ╎    ╎    ╎    ╎    ╎    ╎    7      …sentials.jl:647; ifelse
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   5      …ase/math.jl:911; max
    5╎    ╎    ╎    ╎    ╎    ╎    ╎    5      …sentials.jl:647; ifelse
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   5      …/generic.jl:639; norm
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    5      …/generic.jl:639; norm
    5╎    ╎    ╎    ╎    ╎    ╎    ╎     5      …se/float.jl:610; abs
    1╎    ╎    ╎    ╎    ╎    ╎ 11     …/generic.jl:464; generic_norm2(x::Vecto…
    3╎    ╎    ╎    ╎    ╎    ╎  10     …se/float.jl:635; isinf
    7╎    ╎    ╎    ╎    ╎    ╎   7      …se/float.jl:610; abs
     ╎    ╎    ╎    ╎    ╎    ╎ 6      …/generic.jl:465; generic_norm2(x::Vecto…
     ╎    ╎    ╎    ╎    ╎    ╎  6      …se/array.jl:945; iterate
     ╎    ╎    ╎    ╎    ╎    ╎   6      …se/array.jl:945; iterate
    6╎    ╎    ╎    ╎    ╎    ╎    6      …sentials.jl:13; getindex
   11╎    ╎    ╎    ╎    ╎    ╎ 30     …/generic.jl:467; generic_norm2(x::Vecto…
     ╎    ╎    ╎    ╎    ╎    ╎  13     …se/float.jl:623; isfinite
    2╎    ╎    ╎    ╎    ╎    ╎   2      …se/float.jl:410; -
     ╎    ╎    ╎    ╎    ╎    ╎   11     …se/float.jl:620; isnan
   11╎    ╎    ╎    ╎    ╎    ╎    11     …se/float.jl:535; !=
     ╎    ╎    ╎    ╎    ╎    ╎  6      …perators.jl:587; *
    6╎    ╎    ╎    ╎    ╎    ╎   6      …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎ 4      …/generic.jl:468; generic_norm2(x::Vecto…
     ╎    ╎    ╎    ╎    ╎    ╎  4      …/generic.jl:459; norm_sqr
     ╎    ╎    ╎    ╎    ╎    ╎   4      …e/number.jl:189; abs2
    4╎    ╎    ╎    ╎    ╎    ╎    4      …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎ 17     …/generic.jl:470; generic_norm2(x::Vecto…
    6╎    ╎    ╎    ╎    ╎    ╎  17     …se/array.jl:945; iterate
    8╎    ╎    ╎    ╎    ╎    ╎   8      …sentials.jl:13; getindex
    3╎    ╎    ╎    ╎    ╎    ╎   3      @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎    ╎ 3      …/generic.jl:473; generic_norm2(x::Vecto…
    3╎    ╎    ╎    ╎    ╎    ╎  3      …se/float.jl:409; +
   15╎    ╎    ╎    ╎    ╎  15     …pse_comp.jl:100; eclipse_compute_quantities…
     ╎    ╎    ╎    ╎    ╎  8      …pse_comp.jl:102; eclipse_compute_quantities…
     ╎    ╎    ╎    ╎    ╎   8      …roadcast.jl:911; materialize!
     ╎    ╎    ╎    ╎    ╎    8      …roadcast.jl:914; materialize!
     ╎    ╎    ╎    ╎    ╎     8      …roadcast.jl:956; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎ 8      …roadcast.jl:1003; copyto!
    1╎    ╎    ╎    ╎    ╎    ╎  1      …simdloop.jl:0; macro expansion
    2╎    ╎    ╎    ╎    ╎    ╎  2      …simdloop.jl:75; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎  4      …simdloop.jl:77; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎   4      …roadcast.jl:1004; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎    1      …roadcast.jl:636; getindex
     ╎    ╎    ╎    ╎    ╎    ╎     1      …roadcast.jl:682; _broadcast_getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …roadcast.jl:709; _broadcast_getind…
    1╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎    3      …ensional.jl:698; setindex!
    3╎    ╎    ╎    ╎    ╎    ╎     3      …se/array.jl:1024; setindex!
     ╎    ╎    ╎    ╎    ╎    ╎  1      …simdloop.jl:84; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎   1      …ensional.jl:408; iterate
    1╎    ╎    ╎    ╎    ╎    ╎    1      …ensional.jl:427; __inc
     ╎    ╎    ╎    ╎    ╎  54     …pse_comp.jl:104; eclipse_compute_quantities…
     ╎    ╎    ╎    ╎    ╎   34     …rraymath.jl:24; /(A::Matrix{Float64}, B::F…
     ╎    ╎    ╎    ╎    ╎    34     …roadcast.jl:892; broadcast_preserving_zer…
     ╎    ╎    ╎    ╎    ╎     34     …roadcast.jl:903; materialize
     ╎    ╎    ╎    ╎    ╎    ╎ 34     …roadcast.jl:928; copy
     ╎    ╎    ╎    ╎    ╎    ╎  17     …roadcast.jl:956; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎   17     …roadcast.jl:1003; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎    1      …simdloop.jl:75; macro expansion
    1╎    ╎    ╎    ╎    ╎    ╎     1      @Base/int.jl:83; <
     ╎    ╎    ╎    ╎    ╎    ╎    16     …simdloop.jl:77; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎     16     …roadcast.jl:1004; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 6      …roadcast.jl:636; getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …roadcast.jl:681; _broadcast_getin…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …roadcast.jl:705; _getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …roadcast.jl:675; _broadcast_get…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …ensional.jl:696; getindex
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …sentials.jl:14; getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  5      …roadcast.jl:682; _broadcast_getin…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   5      …roadcast.jl:709; _broadcast_geti…
    5╎    ╎    ╎    ╎    ╎    ╎    ╎    5      …se/float.jl:412; /
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 10     …ensional.jl:698; setindex!
   10╎    ╎    ╎    ╎    ╎    ╎    ╎  10     …se/array.jl:1024; setindex!
     ╎    ╎    ╎    ╎    ╎    ╎  17     …roadcast.jl:223; similar
     ╎    ╎    ╎    ╎    ╎    ╎   17     …roadcast.jl:224; similar
     ╎    ╎    ╎    ╎    ╎    ╎    17     …actarray.jl:876; similar
     ╎    ╎    ╎    ╎    ╎    ╎     17     …actarray.jl:877; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 17     …ase/boot.jl:494; Array
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  17     …ase/boot.jl:487; Array
   17╎    ╎    ╎    ╎    ╎    ╎    ╎   17     …ase/boot.jl:479; Array
     ╎    ╎    ╎    ╎    ╎   20     …roadcast.jl:911; materialize!
     ╎    ╎    ╎    ╎    ╎    20     …roadcast.jl:914; materialize!
     ╎    ╎    ╎    ╎    ╎     20     …roadcast.jl:956; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎ 4      …roadcast.jl:1000; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎  4      …roadcast.jl:983; preprocess
     ╎    ╎    ╎    ╎    ╎    ╎   4      …roadcast.jl:986; preprocess_args
     ╎    ╎    ╎    ╎    ╎    ╎    4      …roadcast.jl:987; preprocess_args
     ╎    ╎    ╎    ╎    ╎    ╎     4      …roadcast.jl:984; preprocess
    1╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …roadcast.jl:977; broadcast_unalias
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 3      …roadcast.jl:676; extrude
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  3      …roadcast.jl:625; newindexer
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   3      …roadcast.jl:626; shapeindexer
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    2      …roadcast.jl:630; _newindexer
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     2      …roadcast.jl:631; _newindexer
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 2      …perators.jl:276; !=
    2╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  2      …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …roadcast.jl:631; _newindexer
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …perators.jl:276; !=
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎ 16     …roadcast.jl:1003; copyto!
    2╎    ╎    ╎    ╎    ╎    ╎  2      …simdloop.jl:0; macro expansion
    5╎    ╎    ╎    ╎    ╎    ╎  5      …simdloop.jl:75; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎  9      …simdloop.jl:77; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎   9      …roadcast.jl:1004; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎    2      …roadcast.jl:636; getindex
     ╎    ╎    ╎    ╎    ╎    ╎     1      …roadcast.jl:681; _broadcast_getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …roadcast.jl:705; _getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …roadcast.jl:675; _broadcast_getin…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …ensional.jl:696; getindex
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …sentials.jl:14; getindex
     ╎    ╎    ╎    ╎    ╎    ╎     1      …roadcast.jl:682; _broadcast_getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …roadcast.jl:709; _broadcast_getind…
    1╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …se/float.jl:409; +
     ╎    ╎    ╎    ╎    ╎    ╎    7      …ensional.jl:698; setindex!
    7╎    ╎    ╎    ╎    ╎    ╎     7      …se/array.jl:1024; setindex!
   23╎    ╎    ╎    ╎    ╎  3407   …pse_comp.jl:109; eclipse_compute_quantities…
    7╎    ╎    ╎    ╎    ╎   7      …actarray.jl:0; getindex
     ╎    ╎    ╎    ╎    ╎   14     …se/array.jl:973; getindex
   10╎    ╎    ╎    ╎    ╎    10     …actarray.jl:700; checkbounds
    4╎    ╎    ╎    ╎    ╎    4      …actarray.jl:702; checkbounds
     ╎    ╎    ╎    ╎    ╎   1931   …se/array.jl:975; getindex
     ╎    ╎    ╎    ╎    ╎    1931   …actarray.jl:831; similar
     ╎    ╎    ╎    ╎    ╎     1931   …se/array.jl:420; similar
     ╎    ╎    ╎    ╎    ╎    ╎ 1931   …ase/boot.jl:486; Array
 1931╎    ╎    ╎    ╎    ╎    ╎  1931   …ase/boot.jl:477; Array
     ╎    ╎    ╎    ╎    ╎   255    …se/array.jl:977; getindex
   46╎    ╎    ╎    ╎    ╎    255    …se/array.jl:368; copyto!
    9╎    ╎    ╎    ╎    ╎     9      …se/array.jl:0; _copyto_impl!(dest::Vecto…
   26╎    ╎    ╎    ╎    ╎     26     …se/array.jl:371; _copyto_impl!(dest::Vec…
     ╎    ╎    ╎    ╎    ╎     4      …se/array.jl:372; _copyto_impl!(dest::Vec…
    4╎    ╎    ╎    ╎    ╎    ╎ 4      …romotion.jl:521; ==
    9╎    ╎    ╎    ╎    ╎     9      …se/array.jl:373; _copyto_impl!(dest::Vec…
     ╎    ╎    ╎    ╎    ╎     42     …se/array.jl:374; _copyto_impl!(dest::Vec…
    5╎    ╎    ╎    ╎    ╎    ╎ 5      …actarray.jl:0; checkbounds
    9╎    ╎    ╎    ╎    ╎    ╎ 9      …actarray.jl:700; checkbounds
     ╎    ╎    ╎    ╎    ╎    ╎ 17     …actarray.jl:702; checkbounds
     ╎    ╎    ╎    ╎    ╎    ╎  17     …actarray.jl:687; checkbounds
     ╎    ╎    ╎    ╎    ╎    ╎   8      …actarray.jl:768; checkindex
     ╎    ╎    ╎    ╎    ╎    ╎    8      …actarray.jl:763; checkindex
    2╎    ╎    ╎    ╎    ╎    ╎     2      @Base/int.jl:86; -
    6╎    ╎    ╎    ╎    ╎    ╎     6      @Base/int.jl:513; <
     ╎    ╎    ╎    ╎    ╎    ╎   9      …actarray.jl:389; eachindex
     ╎    ╎    ╎    ╎    ╎    ╎    9      …actarray.jl:137; axes1
     ╎    ╎    ╎    ╎    ╎    ╎     9      …actarray.jl:98; axes
    9╎    ╎    ╎    ╎    ╎    ╎    ╎ 9      …se/array.jl:191; size
    5╎    ╎    ╎    ╎    ╎    ╎ 5      @Base/int.jl:86; -
     ╎    ╎    ╎    ╎    ╎    ╎ 6      …se/range.jl:5; Colon
     ╎    ╎    ╎    ╎    ╎    ╎  6      …se/range.jl:403; UnitRange
     ╎    ╎    ╎    ╎    ╎    ╎   6      …se/range.jl:414; unitrange_last
     ╎    ╎    ╎    ╎    ╎    ╎    6      …perators.jl:425; >=
    6╎    ╎    ╎    ╎    ╎    ╎     6      @Base/int.jl:514; <=
     ╎    ╎    ╎    ╎    ╎     45     …se/array.jl:375; _copyto_impl!(dest::Vec…
     ╎    ╎    ╎    ╎    ╎    ╎ 29     …actarray.jl:702; checkbounds
     ╎    ╎    ╎    ╎    ╎    ╎  29     …actarray.jl:687; checkbounds
     ╎    ╎    ╎    ╎    ╎    ╎   22     …actarray.jl:768; checkindex
     ╎    ╎    ╎    ╎    ╎    ╎    22     …actarray.jl:763; checkindex
   22╎    ╎    ╎    ╎    ╎    ╎     22     @Base/int.jl:513; <
     ╎    ╎    ╎    ╎    ╎    ╎   7      …actarray.jl:389; eachindex
     ╎    ╎    ╎    ╎    ╎    ╎    7      …actarray.jl:137; axes1
     ╎    ╎    ╎    ╎    ╎    ╎     7      …actarray.jl:98; axes
    7╎    ╎    ╎    ╎    ╎    ╎    ╎ 7      …se/array.jl:191; size
    5╎    ╎    ╎    ╎    ╎    ╎ 5      @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎    ╎ 11     …se/range.jl:5; Colon
     ╎    ╎    ╎    ╎    ╎    ╎  11     …se/range.jl:403; UnitRange
    5╎    ╎    ╎    ╎    ╎    ╎   11     …se/range.jl:414; unitrange_last
     ╎    ╎    ╎    ╎    ╎    ╎    6      …perators.jl:425; >=
    6╎    ╎    ╎    ╎    ╎    ╎     6      @Base/int.jl:514; <=
     ╎    ╎    ╎    ╎    ╎     74     …se/array.jl:376; _copyto_impl!(dest::Vec…
     ╎    ╎    ╎    ╎    ╎    ╎ 8      …se/array.jl:331; unsafe_copyto!
     ╎    ╎    ╎    ╎    ╎    ╎  8      …actarray.jl:1240; pointer
    8╎    ╎    ╎    ╎    ╎    ╎   8      …/pointer.jl:282; +
     ╎    ╎    ╎    ╎    ╎    ╎ 13     …se/array.jl:332; unsafe_copyto!
     ╎    ╎    ╎    ╎    ╎    ╎  13     …actarray.jl:1240; pointer
   13╎    ╎    ╎    ╎    ╎    ╎   13     …/pointer.jl:282; +
     ╎    ╎    ╎    ╎    ╎    ╎ 53     …se/array.jl:337; unsafe_copyto!
   43╎    ╎    ╎    ╎    ╎    ╎  50     …ase/cmem.jl:26; memmove
     ╎    ╎    ╎    ╎    ╎    ╎   7      …sentials.jl:543; cconvert
     ╎    ╎    ╎    ╎    ╎    ╎    7      …e/number.jl:7; convert
     ╎    ╎    ╎    ╎    ╎    ╎     7      …ase/boot.jl:789; UInt64
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 7      …ase/boot.jl:759; toUInt64
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  7      …ase/boot.jl:648; check_top_bit
    7╎    ╎    ╎    ╎    ╎    ╎    ╎   7      …ase/boot.jl:638; is_top_bit_set
    3╎    ╎    ╎    ╎    ╎    ╎  3      @Base/int.jl:88; *
    1╎    ╎    ╎    ╎    ╎   1      …sentials.jl:13; getindex
    8╎    ╎    ╎    ╎    ╎   8      …ial/trig.jl:671; acos(x::Float64)
    3╎    ╎    ╎    ╎    ╎   3      …/generic.jl:595; norm(itr::Vector{Float64}…
   37╎    ╎    ╎    ╎    ╎   37     …geometry.jl:34; calc_proj_dist(p1::Vector{…
   93╎    ╎    ╎    ╎    ╎   1128   …geometry.jl:35; calc_proj_dist(p1::Vector{…
    8╎    ╎    ╎    ╎    ╎    8      …se/float.jl:411; *
   28╎    ╎    ╎    ╎    ╎    28     …ial/trig.jl:0; acos(x::Float64)
    6╎    ╎    ╎    ╎    ╎    6      …ial/trig.jl:671; acos(x::Float64)
     ╎    ╎    ╎    ╎    ╎    3      …ial/trig.jl:695; acos(x::Float64)
    3╎    ╎    ╎    ╎    ╎     3      …se/float.jl:610; abs
     ╎    ╎    ╎    ╎    ╎    7      …ial/trig.jl:696; acos(x::Float64)
     ╎    ╎    ╎    ╎    ╎     7      …perators.jl:425; >=
    7╎    ╎    ╎    ╎    ╎    ╎ 7      …se/float.jl:537; <=
    1╎    ╎    ╎    ╎    ╎    1      …ial/trig.jl:701; acos(x::Float64)
     ╎    ╎    ╎    ╎    ╎    33     …ial/trig.jl:707; acos(x::Float64)
   33╎    ╎    ╎    ╎    ╎     33     …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    152    …ial/trig.jl:708; acos(x::Float64)
     ╎    ╎    ╎    ╎    ╎     152    …ial/trig.jl:391; arc_tRt
     ╎    ╎    ╎    ╎    ╎    ╎ 128    …ial/trig.jl:366; arc_p
    3╎    ╎    ╎    ╎    ╎    ╎  3      …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎  125    …ase/math.jl:186; evalpoly
     ╎    ╎    ╎    ╎    ╎    ╎   125    …ase/math.jl:187; macro expansion
  125╎    ╎    ╎    ╎    ╎    ╎    125    …se/float.jl:414; muladd
     ╎    ╎    ╎    ╎    ╎    ╎ 24     …ial/trig.jl:375; arc_q
     ╎    ╎    ╎    ╎    ╎    ╎  24     …ase/math.jl:186; evalpoly
     ╎    ╎    ╎    ╎    ╎    ╎   24     …ase/math.jl:187; macro expansion
   24╎    ╎    ╎    ╎    ╎    ╎    24     …se/float.jl:414; muladd
   80╎    ╎    ╎    ╎    ╎    80     …ial/trig.jl:709; acos(x::Float64)
     ╎    ╎    ╎    ╎    ╎    5      …ial/trig.jl:720; acos(x::Float64)
     ╎    ╎    ╎    ╎    ╎     5      …ial/trig.jl:668; ACOS_CORRECT_LOWWORD
    5╎    ╎    ╎    ╎    ╎    ╎ 5      …sentials.jl:581; reinterpret
     ╎    ╎    ╎    ╎    ╎    5      …ial/trig.jl:721; acos(x::Float64)
    5╎    ╎    ╎    ╎    ╎     5      …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    72     …ial/trig.jl:722; acos(x::Float64)
   58╎    ╎    ╎    ╎    ╎     58     …se/float.jl:411; *
   14╎    ╎    ╎    ╎    ╎     14     …se/float.jl:409; +
   14╎    ╎    ╎    ╎    ╎    547    …/generic.jl:596; norm
   10╎    ╎    ╎    ╎    ╎     10     …/generic.jl:0; generic_norm2(x::Vector{F…
   18╎    ╎    ╎    ╎    ╎     18     …/generic.jl:595; norm(itr::Vector{Float6…
     ╎    ╎    ╎    ╎    ╎     5      …/generic.jl:596; norm(itr::Vector{Float6…
     ╎    ╎    ╎    ╎    ╎    ╎ 5      …actarray.jl:1220; isempty
    5╎    ╎    ╎    ╎    ╎    ╎  5      …sentials.jl:10; length
   20╎    ╎    ╎    ╎    ╎     20     …/generic.jl:597; norm(itr::Vector{Float6…
     ╎    ╎    ╎    ╎    ╎     480    …/generic.jl:598; norm(itr::Vector{Float6…
    8╎    ╎    ╎    ╎    ╎    ╎ 480    …rc/dense.jl:106; norm2
   22╎    ╎    ╎    ╎    ╎    ╎  22     …se/array.jl:0; generic_norm2(x::Vector…
    5╎    ╎    ╎    ╎    ╎    ╎  5      …e/reduce.jl:0; _mapreduce(f::typeof(Li…
    8╎    ╎    ╎    ╎    ╎    ╎  8      …e/reduce.jl:428; _mapreduce(f::typeof(…
   47╎    ╎    ╎    ╎    ╎    ╎  47     …/generic.jl:0; generic_norm2(x::Vector…
   34╎    ╎    ╎    ╎    ╎    ╎  34     …/generic.jl:462; generic_norm2(x::Vect…
     ╎    ╎    ╎    ╎    ╎    ╎  201    …/generic.jl:463; generic_norm2(x::Vect…
     ╎    ╎    ╎    ╎    ╎    ╎   201    …/generic.jl:527; normInf
     ╎    ╎    ╎    ╎    ╎    ╎    201    …/generic.jl:453; generic_normInf
     ╎    ╎    ╎    ╎    ╎    ╎     201    …educedim.jl:357; mapreduce
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 201    …educedim.jl:357; #mapreduce#821
   14╎    ╎    ╎    ╎    ╎    ╎    ╎  201    …educedim.jl:365; _mapreduce_dim
   11╎    ╎    ╎    ╎    ╎    ╎    ╎   11     …e/reduce.jl:0; _mapreduce(f::typ…
   23╎    ╎    ╎    ╎    ╎    ╎    ╎   23     …e/reduce.jl:428; _mapreduce(f::t…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   5      …e/reduce.jl:429; _mapreduce(f::t…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    5      …/indices.jl:486; LinearIndices
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     5      …actarray.jl:98; axes
    5╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 5      …se/array.jl:191; size
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   7      …e/reduce.jl:436; _mapreduce(f::t…
    7╎    ╎    ╎    ╎    ╎    ╎    ╎    7      @Base/int.jl:83; <
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   50     …e/reduce.jl:440; _mapreduce(f::t…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    11     …ase/math.jl:908; max
   11╎    ╎    ╎    ╎    ╎    ╎    ╎     11     …se/float.jl:410; -
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    17     …ase/math.jl:909; max
   16╎    ╎    ╎    ╎    ╎    ╎    ╎     16     …sentials.jl:647; ifelse
    1╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …oatfuncs.jl:15; signbit
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    10     …ase/math.jl:911; max
   10╎    ╎    ╎    ╎    ╎    ╎    ╎     10     …sentials.jl:647; ifelse
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    12     …/generic.jl:639; norm
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     12     …/generic.jl:639; norm
   12╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 12     …se/float.jl:610; abs
   47╎    ╎    ╎    ╎    ╎    ╎    ╎   47     …e/reduce.jl:441; _mapreduce(f::t…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   8      …e/reduce.jl:442; _mapreduce(f::t…
    8╎    ╎    ╎    ╎    ╎    ╎    ╎    8      …sentials.jl:13; getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   36     …e/reduce.jl:443; _mapreduce(f::t…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    20     …ase/math.jl:909; max
   20╎    ╎    ╎    ╎    ╎    ╎    ╎     20     …sentials.jl:647; ifelse
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    6      …ase/math.jl:911; max
    6╎    ╎    ╎    ╎    ╎    ╎    ╎     6      …sentials.jl:647; ifelse
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    10     …/generic.jl:639; norm
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     10     …/generic.jl:639; norm
   10╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 10     …se/float.jl:610; abs
    4╎    ╎    ╎    ╎    ╎    ╎  38     …/generic.jl:464; generic_norm2(x::Vect…
   19╎    ╎    ╎    ╎    ╎    ╎   34     …se/float.jl:635; isinf
   15╎    ╎    ╎    ╎    ╎    ╎    15     …se/float.jl:610; abs
     ╎    ╎    ╎    ╎    ╎    ╎  16     …/generic.jl:465; generic_norm2(x::Vect…
     ╎    ╎    ╎    ╎    ╎    ╎   16     …se/array.jl:945; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    16     …se/array.jl:945; iterate
   16╎    ╎    ╎    ╎    ╎    ╎     16     …sentials.jl:13; getindex
    7╎    ╎    ╎    ╎    ╎    ╎  26     …/generic.jl:467; generic_norm2(x::Vect…
     ╎    ╎    ╎    ╎    ╎    ╎   7      …se/float.jl:623; isfinite
    4╎    ╎    ╎    ╎    ╎    ╎    4      …se/float.jl:410; -
     ╎    ╎    ╎    ╎    ╎    ╎    3      …se/float.jl:620; isnan
    3╎    ╎    ╎    ╎    ╎    ╎     3      …se/float.jl:535; !=
     ╎    ╎    ╎    ╎    ╎    ╎   12     …perators.jl:587; *
   12╎    ╎    ╎    ╎    ╎    ╎    12     …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎  9      …/generic.jl:468; generic_norm2(x::Vect…
     ╎    ╎    ╎    ╎    ╎    ╎   9      …/generic.jl:459; norm_sqr
     ╎    ╎    ╎    ╎    ╎    ╎    9      …e/number.jl:189; abs2
    9╎    ╎    ╎    ╎    ╎    ╎     9      …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎  49     …/generic.jl:470; generic_norm2(x::Vect…
   19╎    ╎    ╎    ╎    ╎    ╎   49     …se/array.jl:945; iterate
   17╎    ╎    ╎    ╎    ╎    ╎    17     …sentials.jl:13; getindex
   13╎    ╎    ╎    ╎    ╎    ╎    13     @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎    ╎  17     …/generic.jl:473; generic_norm2(x::Vect…
   17╎    ╎    ╎    ╎    ╎    ╎   17     …se/float.jl:409; +
     ╎    ╎    ╎    ╎    ╎    88     …c/matmul.jl:15; dot
     ╎    ╎    ╎    ╎    ╎     6      …src/blas.jl:393; dot
    6╎    ╎    ╎    ╎    ╎    ╎ 6      …sentials.jl:10; length
     ╎    ╎    ╎    ╎    ╎     82     …src/blas.jl:395; dot
   74╎    ╎    ╎    ╎    ╎    ╎ 74     …src/blas.jl:345; dot
     ╎    ╎    ╎    ╎    ╎    ╎ 8      …src/blas.jl:177; vec_pointer_stride
     ╎    ╎    ╎    ╎    ╎    ╎  8      …src/blas.jl:177; vec_pointer_stride
     ╎    ╎    ╎    ╎    ╎    ╎   8      …actarray.jl:1237; pointer
    8╎    ╎    ╎    ╎    ╎    ╎    8      …/pointer.jl:65; unsafe_convert
   12╎    ╎    ╎    ╎    ╎  12     …pse_comp.jl:110; eclipse_compute_quantities…
     ╎    ╎    ╎    ╎    ╎  39     …pse_comp.jl:113; eclipse_compute_quantities…
     ╎    ╎    ╎    ╎    ╎   39     …roadcast.jl:903; materialize
     ╎    ╎    ╎    ╎    ╎    39     …roadcast.jl:928; copy
     ╎    ╎    ╎    ╎    ╎     36     …roadcast.jl:956; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎ 3      …roadcast.jl:1012; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎  3      …se/tuple.jl:482; ==
     ╎    ╎    ╎    ╎    ╎    ╎   3      …se/tuple.jl:486; _eq
     ╎    ╎    ╎    ╎    ╎    ╎    3      …se/range.jl:1134; ==
    3╎    ╎    ╎    ╎    ╎    ╎     3      …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎ 17     …roadcast.jl:1015; copyto!
   17╎    ╎    ╎    ╎    ╎    ╎  17     …ase/boot.jl:477; Array
     ╎    ╎    ╎    ╎    ╎    ╎ 2      …roadcast.jl:1019; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎  2      …terators.jl:1271; partition
     ╎    ╎    ╎    ╎    ╎    ╎   2      …terators.jl:1279; PartitionIterator
     ╎    ╎    ╎    ╎    ╎    ╎    2      …rraymath.jl:41; vec
     ╎    ╎    ╎    ╎    ╎    ╎     2      …pedarray.jl:117; reshape
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 2      …pedarray.jl:112; reshape
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  2      …pedarray.jl:178; _reshape
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   2      …pedarray.jl:193; __reshape
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    2      …se/tuple.jl:291; map
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     2      …inverses.jl:90; SignedMultipli…
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      @Base/int.jl:0; Base.Multiplic…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …inverses.jl:70; Base.Multipli…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …perators.jl:425; >=
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      @Base/int.jl:515; <=
     ╎    ╎    ╎    ╎    ╎    ╎ 8      …roadcast.jl:1021; copyto!
    1╎    ╎    ╎    ╎    ╎    ╎  1      …simdloop.jl:0; macro expansion
    1╎    ╎    ╎    ╎    ╎    ╎  1      …simdloop.jl:72; macro expansion
    1╎    ╎    ╎    ╎    ╎    ╎  2      …simdloop.jl:75; macro expansion
    1╎    ╎    ╎    ╎    ╎    ╎   1      @Base/int.jl:83; <
     ╎    ╎    ╎    ╎    ╎    ╎  4      …simdloop.jl:77; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎   4      …roadcast.jl:1022; macro expansion
    2╎    ╎    ╎    ╎    ╎    ╎    2      …se/array.jl:1021; setindex!
     ╎    ╎    ╎    ╎    ╎    ╎    2      …roadcast.jl:636; getindex
     ╎    ╎    ╎    ╎    ╎    ╎     1      …roadcast.jl:681; _broadcast_getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …roadcast.jl:705; _getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …roadcast.jl:675; _broadcast_getin…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …ensional.jl:696; getindex
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …sentials.jl:14; getindex
     ╎    ╎    ╎    ╎    ╎    ╎     1      …roadcast.jl:682; _broadcast_getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …roadcast.jl:709; _broadcast_getind…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …perators.jl:378; >
    1╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …se/float.jl:536; <
     ╎    ╎    ╎    ╎    ╎    ╎ 1      …roadcast.jl:1025; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎  1      …simdloop.jl:77; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎   1      …roadcast.jl:1026; macro expansion
    1╎    ╎    ╎    ╎    ╎    ╎    1      …se/array.jl:1021; setindex!
     ╎    ╎    ╎    ╎    ╎    ╎ 3      …roadcast.jl:1028; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎  3      …bitarray.jl:356; dumpbitcache
     ╎    ╎    ╎    ╎    ╎    ╎   2      …bitarray.jl:325; copy_to_bitarray_chu…
    2╎    ╎    ╎    ╎    ╎    ╎    2      …se/range.jl:901; iterate
     ╎    ╎    ╎    ╎    ╎    ╎   1      …bitarray.jl:328; copy_to_bitarray_chu…
     ╎    ╎    ╎    ╎    ╎    ╎    1      …se/range.jl:897; iterate
     ╎    ╎    ╎    ╎    ╎    ╎     1      …se/range.jl:672; isempty
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …perators.jl:378; >
    1╎    ╎    ╎    ╎    ╎    ╎    ╎  1      @Base/int.jl:83; <
    2╎    ╎    ╎    ╎    ╎    ╎ 2      @Base/int.jl:0; copyto!
     ╎    ╎    ╎    ╎    ╎     3      …roadcast.jl:223; similar
     ╎    ╎    ╎    ╎    ╎    ╎ 3      …roadcast.jl:226; similar
     ╎    ╎    ╎    ╎    ╎    ╎  3      …actarray.jl:876; similar
     ╎    ╎    ╎    ╎    ╎    ╎   3      …actarray.jl:877; similar
    1╎    ╎    ╎    ╎    ╎    ╎    3      …bitarray.jl:71; BitArray
     ╎    ╎    ╎    ╎    ╎    ╎     2      …bitarray.jl:37; BitMatrix(::UndefIn…
    2╎    ╎    ╎    ╎    ╎    ╎    ╎ 2      …ase/boot.jl:477; Array
   31╎    ╎    ╎    ╎    ╎  108    …pse_comp.jl:115; eclipse_compute_quantities…
     ╎    ╎    ╎    ╎    ╎   1      …se/array.jl:973; getindex
    1╎    ╎    ╎    ╎    ╎    1      …actarray.jl:700; checkbounds
     ╎    ╎    ╎    ╎    ╎   11     …se/array.jl:975; getindex
     ╎    ╎    ╎    ╎    ╎    11     …actarray.jl:831; similar
     ╎    ╎    ╎    ╎    ╎     11     …se/array.jl:420; similar
     ╎    ╎    ╎    ╎    ╎    ╎ 11     …ase/boot.jl:486; Array
   11╎    ╎    ╎    ╎    ╎    ╎  11     …ase/boot.jl:477; Array
    1╎    ╎    ╎    ╎    ╎   1      …roadcast.jl:1343; broadcasted(::typeof(>),…
    1╎    ╎    ╎    ╎    ╎   1      …roadcast.jl:0; materialize(bc::Base.Broadc…
     ╎    ╎    ╎    ╎    ╎   59     …roadcast.jl:903; materialize(bc::Base.Broa…
     ╎    ╎    ╎    ╎    ╎    58     …roadcast.jl:928; copy
     ╎    ╎    ╎    ╎    ╎     54     …roadcast.jl:956; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎ 7      …roadcast.jl:1015; copyto!
    7╎    ╎    ╎    ╎    ╎    ╎  7      …ase/boot.jl:477; Array
     ╎    ╎    ╎    ╎    ╎    ╎ 1      …roadcast.jl:1018; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎  1      …roadcast.jl:983; preprocess
     ╎    ╎    ╎    ╎    ╎    ╎   1      …roadcast.jl:986; preprocess_args
     ╎    ╎    ╎    ╎    ╎    ╎    1      …roadcast.jl:987; preprocess_args
     ╎    ╎    ╎    ╎    ╎    ╎     1      …roadcast.jl:983; preprocess
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …roadcast.jl:986; preprocess_args
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …roadcast.jl:984; preprocess
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …roadcast.jl:977; broadcast_unali…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …actarray.jl:1481; unalias
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …actarray.jl:1516; mightalias
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …actarray.jl:1539; dataids
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …flection.jl:611; objectid
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …flection.jl:617; _objectid
     ╎    ╎    ╎    ╎    ╎    ╎ 1      …roadcast.jl:1019; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎  1      …terators.jl:1271; partition
     ╎    ╎    ╎    ╎    ╎    ╎   1      …terators.jl:1279; PartitionIterator
     ╎    ╎    ╎    ╎    ╎    ╎    1      …rraymath.jl:41; vec
     ╎    ╎    ╎    ╎    ╎    ╎     1      …pedarray.jl:117; reshape
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …pedarray.jl:112; reshape
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …pedarray.jl:178; _reshape
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …pedarray.jl:193; __reshape
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …se/tuple.jl:291; map
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …inverses.jl:90; SignedMultipli…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …inverses.jl:60; Base.Multipli…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      @Base/int.jl:1068; +
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎    ╎ 39     …roadcast.jl:1021; copyto!
    1╎    ╎    ╎    ╎    ╎    ╎  1      @Base/int.jl:0; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎  3      …simdloop.jl:75; macro expansion
    3╎    ╎    ╎    ╎    ╎    ╎   3      @Base/int.jl:83; <
     ╎    ╎    ╎    ╎    ╎    ╎  32     …simdloop.jl:77; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎   29     …roadcast.jl:1022; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎    29     …roadcast.jl:636; getindex
     ╎    ╎    ╎    ╎    ╎    ╎     27     …roadcast.jl:681; _broadcast_getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 27     …roadcast.jl:705; _getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  21     …roadcast.jl:675; _broadcast_getin…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   21     …actarray.jl:1291; getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    21     …actarray.jl:1324; _getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     6      …actarray.jl:1330; _to_linear_i…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 6      …actarray.jl:2957; _sub2ind
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  6      …actarray.jl:98; axes
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   6      …se/tuple.jl:292; map
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    6      …se/range.jl:469; oneto
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     6      …se/range.jl:467; OneTo
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 6      …se/range.jl:454; OneTo
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +2 6      …romotion.jl:532; max
    6╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +3 6      …sentials.jl:647; ifelse
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     15     …bitarray.jl:682; getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 8      …bitarray.jl:674; unsafe_bitge…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  8      …bitarray.jl:126; get_chunks_…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   8      …bitarray.jl:119; _div64
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    8      @Base/int.jl:534; >>
    8╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     8      @Base/int.jl:527; >>
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 7      …bitarray.jl:676; unsafe_bitge…
    2╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  2      …sentials.jl:13; getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  5      …perators.jl:276; !=
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   5      @Base/int.jl:518; ==
    5╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    5      …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  6      …roadcast.jl:706; _getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   3      …roadcast.jl:681; _broadcast_geti…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    3      …roadcast.jl:705; _getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     3      …roadcast.jl:675; _broadcast_ge…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 3      …ensional.jl:696; getindex
    3╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  3      …sentials.jl:14; getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   3      …roadcast.jl:682; _broadcast_geti…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    3      …roadcast.jl:709; _broadcast_get…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     3      …perators.jl:378; >
    3╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 3      …se/float.jl:536; <
     ╎    ╎    ╎    ╎    ╎    ╎     2      …roadcast.jl:682; _broadcast_getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 2      …roadcast.jl:709; _broadcast_getind…
    2╎    ╎    ╎    ╎    ╎    ╎    ╎  2      …ase/bool.jl:38; &
     ╎    ╎    ╎    ╎    ╎    ╎   3      …roadcast.jl:1023; macro expansion
    3╎    ╎    ╎    ╎    ╎    ╎    3      @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎    ╎  3      …simdloop.jl:84; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎   3      …enerator.jl:47; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    3      none:0; #30
     ╎    ╎    ╎    ╎    ╎    ╎     3      …ensional.jl:667; skip_len_I
    3╎    ╎    ╎    ╎    ╎    ╎    ╎ 3      @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎    ╎ 6      …roadcast.jl:1028; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎  6      …bitarray.jl:356; dumpbitcache
     ╎    ╎    ╎    ╎    ╎    ╎   1      …bitarray.jl:285; copy_to_bitarray_chu…
    1╎    ╎    ╎    ╎    ╎    ╎    1      @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎    ╎   4      …bitarray.jl:320; copy_to_bitarray_chu…
     ╎    ╎    ╎    ╎    ╎    ╎    3      …bitarray.jl:277; pack8bools
     ╎    ╎    ╎    ╎    ╎    ╎     1      @Base/int.jl:538; >>>
    1╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      @Base/int.jl:530; >>>
    2╎    ╎    ╎    ╎    ╎    ╎     2      @Base/int.jl:372; |
     ╎    ╎    ╎    ╎    ╎    ╎    1      @Base/int.jl:536; <<
    1╎    ╎    ╎    ╎    ╎    ╎     1      @Base/int.jl:529; <<
     ╎    ╎    ╎    ╎    ╎    ╎   1      …bitarray.jl:324; copy_to_bitarray_chu…
    1╎    ╎    ╎    ╎    ╎    ╎    1      @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎     4      …roadcast.jl:223; similar
     ╎    ╎    ╎    ╎    ╎    ╎ 4      …roadcast.jl:226; similar
     ╎    ╎    ╎    ╎    ╎    ╎  4      …actarray.jl:876; similar
     ╎    ╎    ╎    ╎    ╎    ╎   4      …actarray.jl:877; similar
     ╎    ╎    ╎    ╎    ╎    ╎    4      …bitarray.jl:71; BitArray
     ╎    ╎    ╎    ╎    ╎    ╎     4      …bitarray.jl:37; BitMatrix(::UndefIn…
    4╎    ╎    ╎    ╎    ╎    ╎    ╎ 4      …ase/boot.jl:477; Array
    1╎    ╎    ╎    ╎    ╎    1      …simdloop.jl:0; copy
    1╎    ╎    ╎    ╎    ╎   1      …se/float.jl:412; /(x::Float64, y::Float64)
    1╎    ╎    ╎    ╎    ╎   1      …ial/trig.jl:502; atan(x::Float64)
     ╎    ╎    ╎    ╎    ╎   1      …ial/trig.jl:522; atan(x::Float64)
    1╎    ╎    ╎    ╎    ╎    1      …se/float.jl:536; <
     ╎    ╎    ╎    ╎    ╎   1      …ial/trig.jl:527; atan(x::Float64)
     ╎    ╎    ╎    ╎    ╎    1      …ial/trig.jl:499; atan_pq
     ╎    ╎    ╎    ╎    ╎     1      …ial/trig.jl:480; atan_p
     ╎    ╎    ╎    ╎    ╎    ╎ 1      …ase/math.jl:186; evalpoly
     ╎    ╎    ╎    ╎    ╎    ╎  1      …ase/math.jl:187; macro expansion
    1╎    ╎    ╎    ╎    ╎    ╎   1      …se/float.jl:414; muladd
    1╎    ╎    ╎    ╎    ╎  50     …pse_comp.jl:118; eclipse_compute_quantities…
    4╎    ╎    ╎    ╎    ╎   4      …se/array.jl:1024; setindex!
     ╎    ╎    ╎    ╎    ╎   6      …subarray.jl:183; view
     ╎    ╎    ╎    ╎    ╎    1      …/indices.jl:345; to_indices
     ╎    ╎    ╎    ╎    ╎     1      …ensional.jl:869; to_indices
     ╎    ╎    ╎    ╎    ╎    ╎ 1      …ensional.jl:872; _maybe_linear_logical_…
     ╎    ╎    ╎    ╎    ╎    ╎  1      …ensional.jl:789; LogicalIndex
     ╎    ╎    ╎    ╎    ╎    ╎   1      …ensional.jl:785; LogicalIndex
     ╎    ╎    ╎    ╎    ╎    ╎    1      …educedim.jl:439; count
     ╎    ╎    ╎    ╎    ╎    ╎     1      …educedim.jl:439; #count#824
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …educedim.jl:440; count
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …educedim.jl:440; #count#825
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …bitarray.jl:1454; _count
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …bitarray.jl:1446; bitcount
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …bitarray.jl:1449; #bitcount#353
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      @Base/int.jl:415; count_ones
     ╎    ╎    ╎    ╎    ╎    5      …se/tuple.jl:291; map
     ╎    ╎    ╎    ╎    ╎     5      …subarray.jl:183; #214
     ╎    ╎    ╎    ╎    ╎    ╎ 5      …actarray.jl:1481; unalias
     ╎    ╎    ╎    ╎    ╎    ╎  5      …actarray.jl:1516; mightalias
     ╎    ╎    ╎    ╎    ╎    ╎   5      …actarray.jl:1539; dataids
     ╎    ╎    ╎    ╎    ╎    ╎    5      …flection.jl:611; objectid
    5╎    ╎    ╎    ╎    ╎    ╎     5      …flection.jl:617; _objectid
     ╎    ╎    ╎    ╎    ╎   24     …subarray.jl:186; view
     ╎    ╎    ╎    ╎    ╎    2      …subarray.jl:126; _maybe_reshape_parent
     ╎    ╎    ╎    ╎    ╎     2      …pedarray.jl:142; reshape
     ╎    ╎    ╎    ╎    ╎    ╎ 2      …pedarray.jl:111; reshape
    2╎    ╎    ╎    ╎    ╎    ╎  2      …pedarray.jl:51; reshape
     ╎    ╎    ╎    ╎    ╎    22     …subarray.jl:223; unsafe_view
     ╎    ╎    ╎    ╎    ╎     22     …subarray.jl:28; SubArray
     ╎    ╎    ╎    ╎    ╎    ╎ 22     …ensional.jl:855; ensure_indexable
     ╎    ╎    ╎    ╎    ╎    ╎  22     …ensional.jl:792; collect
     ╎    ╎    ╎    ╎    ╎    ╎   7      …se/array.jl:839; collect(itr::Base.Ge…
     ╎    ╎    ╎    ╎    ╎    ╎    7      …se/array.jl:723; _array_for
     ╎    ╎    ╎    ╎    ╎    ╎     7      …actarray.jl:876; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 7      …actarray.jl:877; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  7      …ase/boot.jl:486; Array
    7╎    ╎    ╎    ╎    ╎    ╎    ╎   7      …ase/boot.jl:477; Array
     ╎    ╎    ╎    ╎    ╎    ╎   15     …se/array.jl:844; collect(itr::Base.Ge…
     ╎    ╎    ╎    ╎    ╎    ╎    15     …se/array.jl:870; collect_to_with_fir…
     ╎    ╎    ╎    ╎    ╎    ╎     9      …se/array.jl:892; collect_to!(dest::…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 9      …enerator.jl:44; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …ensional.jl:836; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      @Base/int.jl:518; ==
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  2      …ensional.jl:841; iterate
    2╎    ╎    ╎    ╎    ╎    ╎    ╎   2      @Base/int.jl:441; trailing_zeros
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …ensional.jl:842; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …bitarray.jl:121; _blsr
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    1      @Base/int.jl:347; &
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  5      …ensional.jl:843; iterate
    5╎    ╎    ╎    ╎    ╎    ╎    ╎   5      @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎    ╎     1      …se/array.jl:896; collect_to!(dest::…
    1╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …se/array.jl:1021; setindex!
     ╎    ╎    ╎    ╎    ╎    ╎     5      …se/array.jl:897; collect_to!(dest::…
    5╎    ╎    ╎    ╎    ╎    ╎    ╎ 5      @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎   15     …atistics.jl:174; mean
     ╎    ╎    ╎    ╎    ╎    15     …atistics.jl:174; #mean#2
     ╎    ╎    ╎    ╎    ╎     15     …atistics.jl:187; _mean(f::typeof(identit…
     ╎    ╎    ╎    ╎    ╎    ╎ 15     …educedim.jl:1011; sum
     ╎    ╎    ╎    ╎    ╎    ╎  15     …educedim.jl:1011; #sum#829
     ╎    ╎    ╎    ╎    ╎    ╎   15     …educedim.jl:1015; _sum
     ╎    ╎    ╎    ╎    ╎    ╎    15     …educedim.jl:1015; #_sum#831
     ╎    ╎    ╎    ╎    ╎    ╎     15     …educedim.jl:357; mapreduce
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 15     …educedim.jl:357; #mapreduce#821
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  15     …educedim.jl:365; _mapreduce_dim
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   15     …e/reduce.jl:453; _mapreduce
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    15     …e/reduce.jl:175; mapfoldl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     15     …e/reduce.jl:175; #mapfoldl#298
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 15     …e/reduce.jl:44; mapfoldl_impl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  15     …e/reduce.jl:48; foldl_impl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   15     …e/reduce.jl:60; _foldl_impl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    11     …actarray.jl:1215; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     11     …se/range.jl:901; iterate
   11╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 11     …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    4      …actarray.jl:1217; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     4      …subarray.jl:290; getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 4      …subarray.jl:268; reindex
    4╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +2 4      …sentials.jl:13; getindex
     ╎    ╎    ╎    ╎    ╎  197    …pse_comp.jl:121; eclipse_compute_quantities…
    2╎    ╎    ╎    ╎    ╎   2      …se/array.jl:1024; setindex!
     ╎    ╎    ╎    ╎    ╎   130    …roadcast.jl:903; materialize
     ╎    ╎    ╎    ╎    ╎    130    …roadcast.jl:928; copy
     ╎    ╎    ╎    ╎    ╎     97     …roadcast.jl:956; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎ 4      …roadcast.jl:1000; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎  4      …roadcast.jl:983; preprocess
     ╎    ╎    ╎    ╎    ╎    ╎   4      …roadcast.jl:986; preprocess_args
     ╎    ╎    ╎    ╎    ╎    ╎    4      …roadcast.jl:984; preprocess
     ╎    ╎    ╎    ╎    ╎    ╎     4      …roadcast.jl:977; broadcast_unalias
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 4      …actarray.jl:1481; unalias
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  4      …actarray.jl:1516; mightalias
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   3      …actarray.jl:1523; _isdisjoint
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    3      …perators.jl:276; !=
    3╎    ╎    ╎    ╎    ╎    ╎    ╎     3      …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …actarray.jl:1540; dataids
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …actarray.jl:1237; pointer
    1╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …/pointer.jl:65; unsafe_convert
     ╎    ╎    ╎    ╎    ╎    ╎ 93     …roadcast.jl:1003; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎  91     …simdloop.jl:77; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎   91     …roadcast.jl:1004; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎    91     …roadcast.jl:636; getindex
     ╎    ╎    ╎    ╎    ╎    ╎     4      …roadcast.jl:681; _broadcast_getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 4      …roadcast.jl:705; _getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  4      …roadcast.jl:675; _broadcast_getin…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   4      …ensional.jl:696; getindex
    4╎    ╎    ╎    ╎    ╎    ╎    ╎    4      …sentials.jl:14; getindex
     ╎    ╎    ╎    ╎    ╎    ╎     87     …roadcast.jl:682; _broadcast_getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 87     …roadcast.jl:709; _broadcast_getind…
   87╎    ╎    ╎    ╎    ╎    ╎    ╎  87     …sentials.jl:13; getindex
     ╎    ╎    ╎    ╎    ╎    ╎  2      …simdloop.jl:78; macro expansion
    2╎    ╎    ╎    ╎    ╎    ╎   2      @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎     33     …roadcast.jl:223; similar
     ╎    ╎    ╎    ╎    ╎    ╎ 33     …roadcast.jl:224; similar
     ╎    ╎    ╎    ╎    ╎    ╎  33     …actarray.jl:876; similar
     ╎    ╎    ╎    ╎    ╎    ╎   33     …actarray.jl:877; similar
     ╎    ╎    ╎    ╎    ╎    ╎    33     …ase/boot.jl:494; Array
     ╎    ╎    ╎    ╎    ╎    ╎     33     …ase/boot.jl:487; Array
   33╎    ╎    ╎    ╎    ╎    ╎    ╎ 33     …ase/boot.jl:479; Array
     ╎    ╎    ╎    ╎    ╎   1      …subarray.jl:183; view
     ╎    ╎    ╎    ╎    ╎    1      …se/tuple.jl:291; map
     ╎    ╎    ╎    ╎    ╎     1      …subarray.jl:183; #214
     ╎    ╎    ╎    ╎    ╎    ╎ 1      …actarray.jl:1481; unalias
     ╎    ╎    ╎    ╎    ╎    ╎  1      …actarray.jl:1516; mightalias
     ╎    ╎    ╎    ╎    ╎    ╎   1      …actarray.jl:1539; dataids
     ╎    ╎    ╎    ╎    ╎    ╎    1      …flection.jl:611; objectid
    1╎    ╎    ╎    ╎    ╎    ╎     1      …flection.jl:617; _objectid
     ╎    ╎    ╎    ╎    ╎   36     …subarray.jl:186; view
     ╎    ╎    ╎    ╎    ╎    13     …subarray.jl:126; _maybe_reshape_parent
     ╎    ╎    ╎    ╎    ╎     13     …pedarray.jl:142; reshape
     ╎    ╎    ╎    ╎    ╎    ╎ 13     …pedarray.jl:111; reshape
   13╎    ╎    ╎    ╎    ╎    ╎  13     …pedarray.jl:51; reshape
     ╎    ╎    ╎    ╎    ╎    23     …subarray.jl:223; unsafe_view
     ╎    ╎    ╎    ╎    ╎     23     …subarray.jl:28; SubArray
     ╎    ╎    ╎    ╎    ╎    ╎ 23     …ensional.jl:855; ensure_indexable
     ╎    ╎    ╎    ╎    ╎    ╎  23     …ensional.jl:792; collect
     ╎    ╎    ╎    ╎    ╎    ╎   9      …se/array.jl:839; collect(itr::Base.Ge…
     ╎    ╎    ╎    ╎    ╎    ╎    9      …se/array.jl:723; _array_for
     ╎    ╎    ╎    ╎    ╎    ╎     9      …actarray.jl:876; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 9      …actarray.jl:877; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  9      …ase/boot.jl:486; Array
    9╎    ╎    ╎    ╎    ╎    ╎    ╎   9      …ase/boot.jl:477; Array
     ╎    ╎    ╎    ╎    ╎    ╎   14     …se/array.jl:844; collect(itr::Base.Ge…
     ╎    ╎    ╎    ╎    ╎    ╎    14     …se/array.jl:870; collect_to_with_fir…
     ╎    ╎    ╎    ╎    ╎    ╎     8      …se/array.jl:892; collect_to!(dest::…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 8      …enerator.jl:44; iterate
    1╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …ensional.jl:836; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  2      …ensional.jl:841; iterate
    2╎    ╎    ╎    ╎    ╎    ╎    ╎   2      @Base/int.jl:441; trailing_zeros
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …ensional.jl:842; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …bitarray.jl:121; _blsr
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    1      @Base/int.jl:347; &
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  4      …ensional.jl:843; iterate
    4╎    ╎    ╎    ╎    ╎    ╎    ╎   4      @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎    ╎     2      …se/array.jl:896; collect_to!(dest::…
    2╎    ╎    ╎    ╎    ╎    ╎    ╎ 2      …se/array.jl:1021; setindex!
     ╎    ╎    ╎    ╎    ╎    ╎     4      …se/array.jl:897; collect_to!(dest::…
    4╎    ╎    ╎    ╎    ╎    ╎    ╎ 4      @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎   28     …atistics.jl:174; mean
     ╎    ╎    ╎    ╎    ╎    28     …atistics.jl:174; #mean#2
     ╎    ╎    ╎    ╎    ╎     28     …atistics.jl:187; _mean(f::typeof(identit…
     ╎    ╎    ╎    ╎    ╎    ╎ 28     …educedim.jl:1011; sum
     ╎    ╎    ╎    ╎    ╎    ╎  28     …educedim.jl:1011; #sum#829
     ╎    ╎    ╎    ╎    ╎    ╎   28     …educedim.jl:1015; _sum
     ╎    ╎    ╎    ╎    ╎    ╎    28     …educedim.jl:1015; #_sum#831
     ╎    ╎    ╎    ╎    ╎    ╎     28     …educedim.jl:357; mapreduce
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 28     …educedim.jl:357; #mapreduce#821
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  28     …educedim.jl:365; _mapreduce_dim
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   28     …e/reduce.jl:453; _mapreduce
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    28     …e/reduce.jl:175; mapfoldl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     28     …e/reduce.jl:175; #mapfoldl#298
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 28     …e/reduce.jl:44; mapfoldl_impl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  28     …e/reduce.jl:48; foldl_impl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   28     …e/reduce.jl:60; _foldl_impl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    18     …actarray.jl:1215; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     18     …se/range.jl:901; iterate
   18╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 18     …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    10     …actarray.jl:1217; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     10     …subarray.jl:290; getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 10     …subarray.jl:268; reindex
   10╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +2 10     …sentials.jl:13; getindex
     ╎    ╎    ╎    ╎    ╎  153    …pse_comp.jl:122; eclipse_compute_quantities…
     ╎    ╎    ╎    ╎    ╎   70     …roadcast.jl:903; materialize
     ╎    ╎    ╎    ╎    ╎    70     …roadcast.jl:928; copy
     ╎    ╎    ╎    ╎    ╎     36     …roadcast.jl:956; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎ 36     …roadcast.jl:1003; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎  2      …simdloop.jl:75; macro expansion
    2╎    ╎    ╎    ╎    ╎    ╎   2      @Base/int.jl:83; <
     ╎    ╎    ╎    ╎    ╎    ╎  32     …simdloop.jl:77; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎   32     …roadcast.jl:1004; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎    29     …roadcast.jl:636; getindex
     ╎    ╎    ╎    ╎    ╎    ╎     4      …roadcast.jl:681; _broadcast_getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 4      …roadcast.jl:705; _getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  4      …roadcast.jl:675; _broadcast_getin…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   4      …ensional.jl:696; getindex
    4╎    ╎    ╎    ╎    ╎    ╎    ╎    4      …sentials.jl:14; getindex
     ╎    ╎    ╎    ╎    ╎    ╎     25     …roadcast.jl:682; _broadcast_getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 25     …roadcast.jl:709; _broadcast_getind…
   25╎    ╎    ╎    ╎    ╎    ╎    ╎  25     …sentials.jl:13; getindex
     ╎    ╎    ╎    ╎    ╎    ╎    3      …ensional.jl:698; setindex!
    3╎    ╎    ╎    ╎    ╎    ╎     3      …se/array.jl:1024; setindex!
     ╎    ╎    ╎    ╎    ╎    ╎  2      …simdloop.jl:78; macro expansion
    2╎    ╎    ╎    ╎    ╎    ╎   2      @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎     34     …roadcast.jl:223; similar
     ╎    ╎    ╎    ╎    ╎    ╎ 34     …roadcast.jl:224; similar
     ╎    ╎    ╎    ╎    ╎    ╎  34     …actarray.jl:876; similar
     ╎    ╎    ╎    ╎    ╎    ╎   34     …actarray.jl:877; similar
     ╎    ╎    ╎    ╎    ╎    ╎    34     …ase/boot.jl:494; Array
     ╎    ╎    ╎    ╎    ╎    ╎     34     …ase/boot.jl:487; Array
   34╎    ╎    ╎    ╎    ╎    ╎    ╎ 34     …ase/boot.jl:479; Array
     ╎    ╎    ╎    ╎    ╎   5      …subarray.jl:183; view
     ╎    ╎    ╎    ╎    ╎    2      …/indices.jl:345; to_indices
     ╎    ╎    ╎    ╎    ╎     2      …ensional.jl:869; to_indices
     ╎    ╎    ╎    ╎    ╎    ╎ 2      …ensional.jl:872; _maybe_linear_logical_…
     ╎    ╎    ╎    ╎    ╎    ╎  2      …ensional.jl:789; LogicalIndex
     ╎    ╎    ╎    ╎    ╎    ╎   2      …ensional.jl:785; LogicalIndex
     ╎    ╎    ╎    ╎    ╎    ╎    2      …educedim.jl:439; count
     ╎    ╎    ╎    ╎    ╎    ╎     2      …educedim.jl:439; #count#824
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 2      …educedim.jl:440; count
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  2      …educedim.jl:440; #count#825
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   2      …bitarray.jl:1454; _count
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    2      …bitarray.jl:1446; bitcount
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …bitarray.jl:1448; #bitcount#353
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …se/range.jl:897; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …se/range.jl:672; isempty
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …perators.jl:378; >
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      @Base/int.jl:83; <
    1╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …bitarray.jl:1450; #bitcount#353
     ╎    ╎    ╎    ╎    ╎    3      …se/tuple.jl:291; map
     ╎    ╎    ╎    ╎    ╎     3      …subarray.jl:183; #214
     ╎    ╎    ╎    ╎    ╎    ╎ 3      …actarray.jl:1481; unalias
     ╎    ╎    ╎    ╎    ╎    ╎  3      …actarray.jl:1516; mightalias
     ╎    ╎    ╎    ╎    ╎    ╎   3      …actarray.jl:1539; dataids
     ╎    ╎    ╎    ╎    ╎    ╎    3      …flection.jl:611; objectid
    3╎    ╎    ╎    ╎    ╎    ╎     3      …flection.jl:617; _objectid
     ╎    ╎    ╎    ╎    ╎   60     …subarray.jl:186; view
     ╎    ╎    ╎    ╎    ╎    35     …subarray.jl:126; _maybe_reshape_parent
     ╎    ╎    ╎    ╎    ╎     35     …pedarray.jl:142; reshape
     ╎    ╎    ╎    ╎    ╎    ╎ 35     …pedarray.jl:111; reshape
   35╎    ╎    ╎    ╎    ╎    ╎  35     …pedarray.jl:51; reshape
     ╎    ╎    ╎    ╎    ╎    25     …subarray.jl:223; unsafe_view
     ╎    ╎    ╎    ╎    ╎     25     …subarray.jl:28; SubArray
     ╎    ╎    ╎    ╎    ╎    ╎ 25     …ensional.jl:855; ensure_indexable
     ╎    ╎    ╎    ╎    ╎    ╎  25     …ensional.jl:792; collect
     ╎    ╎    ╎    ╎    ╎    ╎   6      …se/array.jl:839; collect(itr::Base.Ge…
     ╎    ╎    ╎    ╎    ╎    ╎    6      …se/array.jl:723; _array_for
     ╎    ╎    ╎    ╎    ╎    ╎     6      …actarray.jl:876; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 6      …actarray.jl:877; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  6      …ase/boot.jl:486; Array
    6╎    ╎    ╎    ╎    ╎    ╎    ╎   6      …ase/boot.jl:477; Array
     ╎    ╎    ╎    ╎    ╎    ╎   19     …se/array.jl:844; collect(itr::Base.Ge…
     ╎    ╎    ╎    ╎    ╎    ╎    19     …se/array.jl:870; collect_to_with_fir…
     ╎    ╎    ╎    ╎    ╎    ╎     10     …se/array.jl:892; collect_to!(dest::…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 10     …enerator.jl:44; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …ensional.jl:836; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      @Base/int.jl:518; ==
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  3      …ensional.jl:842; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   3      …bitarray.jl:121; _blsr
    3╎    ╎    ╎    ╎    ╎    ╎    ╎    3      @Base/int.jl:347; &
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  6      …ensional.jl:843; iterate
    6╎    ╎    ╎    ╎    ╎    ╎    ╎   6      @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎    ╎     3      …se/array.jl:896; collect_to!(dest::…
    3╎    ╎    ╎    ╎    ╎    ╎    ╎ 3      …se/array.jl:1021; setindex!
     ╎    ╎    ╎    ╎    ╎    ╎     5      …se/array.jl:897; collect_to!(dest::…
    5╎    ╎    ╎    ╎    ╎    ╎    ╎ 5      @Base/int.jl:87; +
    1╎    ╎    ╎    ╎    ╎    ╎     1      …se/array.jl:903; collect_to!(dest::…
     ╎    ╎    ╎    ╎    ╎   18     …atistics.jl:174; mean
     ╎    ╎    ╎    ╎    ╎    18     …atistics.jl:174; #mean#2
     ╎    ╎    ╎    ╎    ╎     18     …atistics.jl:187; _mean(f::typeof(identit…
     ╎    ╎    ╎    ╎    ╎    ╎ 18     …educedim.jl:1011; sum
     ╎    ╎    ╎    ╎    ╎    ╎  18     …educedim.jl:1011; #sum#829
     ╎    ╎    ╎    ╎    ╎    ╎   18     …educedim.jl:1015; _sum
     ╎    ╎    ╎    ╎    ╎    ╎    18     …educedim.jl:1015; #_sum#831
     ╎    ╎    ╎    ╎    ╎    ╎     18     …educedim.jl:357; mapreduce
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 18     …educedim.jl:357; #mapreduce#821
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  18     …educedim.jl:365; _mapreduce_dim
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   18     …e/reduce.jl:453; _mapreduce
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    18     …e/reduce.jl:175; mapfoldl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     18     …e/reduce.jl:175; #mapfoldl#298
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 18     …e/reduce.jl:44; mapfoldl_impl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  18     …e/reduce.jl:48; foldl_impl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   18     …e/reduce.jl:60; _foldl_impl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    14     …actarray.jl:1215; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     14     …se/range.jl:901; iterate
   14╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 14     …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    4      …actarray.jl:1217; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     4      …subarray.jl:290; getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 4      …subarray.jl:268; reindex
    4╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +2 4      …sentials.jl:13; getindex
     ╎    ╎    ╎    ╎    ╎  169    …pse_comp.jl:123; eclipse_compute_quantities…
    1╎    ╎    ╎    ╎    ╎   1      …se/array.jl:0; setindex!
     ╎    ╎    ╎    ╎    ╎   90     …roadcast.jl:903; materialize
     ╎    ╎    ╎    ╎    ╎    90     …roadcast.jl:928; copy
     ╎    ╎    ╎    ╎    ╎     40     …roadcast.jl:956; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎ 39     …roadcast.jl:1003; copyto!
    1╎    ╎    ╎    ╎    ╎    ╎  1      …simdloop.jl:0; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎  2      …simdloop.jl:75; macro expansion
    2╎    ╎    ╎    ╎    ╎    ╎   2      @Base/int.jl:83; <
     ╎    ╎    ╎    ╎    ╎    ╎  30     …simdloop.jl:77; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎   30     …roadcast.jl:1004; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎    25     …roadcast.jl:636; getindex
     ╎    ╎    ╎    ╎    ╎    ╎     5      …roadcast.jl:681; _broadcast_getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 5      …roadcast.jl:705; _getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  5      …roadcast.jl:675; _broadcast_getin…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   5      …ensional.jl:696; getindex
    5╎    ╎    ╎    ╎    ╎    ╎    ╎    5      …sentials.jl:14; getindex
     ╎    ╎    ╎    ╎    ╎    ╎     20     …roadcast.jl:682; _broadcast_getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 20     …roadcast.jl:709; _broadcast_getind…
   20╎    ╎    ╎    ╎    ╎    ╎    ╎  20     …sentials.jl:13; getindex
     ╎    ╎    ╎    ╎    ╎    ╎    5      …ensional.jl:698; setindex!
    5╎    ╎    ╎    ╎    ╎    ╎     5      …se/array.jl:1024; setindex!
     ╎    ╎    ╎    ╎    ╎    ╎  6      …simdloop.jl:78; macro expansion
    6╎    ╎    ╎    ╎    ╎    ╎   6      @Base/int.jl:87; +
    1╎    ╎    ╎    ╎    ╎    ╎ 1      …simdloop.jl:0; copyto!
     ╎    ╎    ╎    ╎    ╎     50     …roadcast.jl:223; similar
     ╎    ╎    ╎    ╎    ╎    ╎ 50     …roadcast.jl:224; similar
     ╎    ╎    ╎    ╎    ╎    ╎  50     …actarray.jl:876; similar
     ╎    ╎    ╎    ╎    ╎    ╎   50     …actarray.jl:877; similar
     ╎    ╎    ╎    ╎    ╎    ╎    50     …ase/boot.jl:494; Array
     ╎    ╎    ╎    ╎    ╎    ╎     50     …ase/boot.jl:487; Array
   50╎    ╎    ╎    ╎    ╎    ╎    ╎ 50     …ase/boot.jl:479; Array
     ╎    ╎    ╎    ╎    ╎   2      …subarray.jl:183; view
     ╎    ╎    ╎    ╎    ╎    1      …/indices.jl:345; to_indices
     ╎    ╎    ╎    ╎    ╎     1      …ensional.jl:869; to_indices
     ╎    ╎    ╎    ╎    ╎    ╎ 1      …ensional.jl:872; _maybe_linear_logical_…
     ╎    ╎    ╎    ╎    ╎    ╎  1      …ensional.jl:789; LogicalIndex
     ╎    ╎    ╎    ╎    ╎    ╎   1      …ensional.jl:785; LogicalIndex
     ╎    ╎    ╎    ╎    ╎    ╎    1      …educedim.jl:439; count
     ╎    ╎    ╎    ╎    ╎    ╎     1      …educedim.jl:439; #count#824
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …educedim.jl:440; count
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …educedim.jl:440; #count#825
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …bitarray.jl:1454; _count
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …bitarray.jl:1446; bitcount
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …bitarray.jl:1449; #bitcount#353
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      @Base/int.jl:415; count_ones
     ╎    ╎    ╎    ╎    ╎    1      …se/tuple.jl:291; map
     ╎    ╎    ╎    ╎    ╎     1      …subarray.jl:183; #214
     ╎    ╎    ╎    ╎    ╎    ╎ 1      …actarray.jl:1481; unalias
     ╎    ╎    ╎    ╎    ╎    ╎  1      …actarray.jl:1516; mightalias
     ╎    ╎    ╎    ╎    ╎    ╎   1      …actarray.jl:1539; dataids
     ╎    ╎    ╎    ╎    ╎    ╎    1      …flection.jl:611; objectid
    1╎    ╎    ╎    ╎    ╎    ╎     1      …flection.jl:617; _objectid
     ╎    ╎    ╎    ╎    ╎   57     …subarray.jl:186; view
     ╎    ╎    ╎    ╎    ╎    35     …subarray.jl:126; _maybe_reshape_parent
     ╎    ╎    ╎    ╎    ╎     35     …pedarray.jl:142; reshape
     ╎    ╎    ╎    ╎    ╎    ╎ 35     …pedarray.jl:111; reshape
   35╎    ╎    ╎    ╎    ╎    ╎  35     …pedarray.jl:51; reshape
     ╎    ╎    ╎    ╎    ╎    22     …subarray.jl:223; unsafe_view
     ╎    ╎    ╎    ╎    ╎     22     …subarray.jl:28; SubArray
     ╎    ╎    ╎    ╎    ╎    ╎ 22     …ensional.jl:855; ensure_indexable
     ╎    ╎    ╎    ╎    ╎    ╎  22     …ensional.jl:792; collect
     ╎    ╎    ╎    ╎    ╎    ╎   8      …se/array.jl:839; collect(itr::Base.Ge…
     ╎    ╎    ╎    ╎    ╎    ╎    8      …se/array.jl:723; _array_for
     ╎    ╎    ╎    ╎    ╎    ╎     8      …actarray.jl:876; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 8      …actarray.jl:877; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  8      …ase/boot.jl:486; Array
    8╎    ╎    ╎    ╎    ╎    ╎    ╎   8      …ase/boot.jl:477; Array
     ╎    ╎    ╎    ╎    ╎    ╎   14     …se/array.jl:844; collect(itr::Base.Ge…
     ╎    ╎    ╎    ╎    ╎    ╎    14     …se/array.jl:870; collect_to_with_fir…
     ╎    ╎    ╎    ╎    ╎    ╎     9      …se/array.jl:892; collect_to!(dest::…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 9      …enerator.jl:44; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …ensional.jl:836; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      @Base/int.jl:518; ==
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  4      …ensional.jl:841; iterate
    4╎    ╎    ╎    ╎    ╎    ╎    ╎   4      @Base/int.jl:441; trailing_zeros
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …ensional.jl:842; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …bitarray.jl:121; _blsr
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    1      @Base/int.jl:347; &
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  3      …ensional.jl:843; iterate
    3╎    ╎    ╎    ╎    ╎    ╎    ╎   3      @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎    ╎     5      …se/array.jl:897; collect_to!(dest::…
    5╎    ╎    ╎    ╎    ╎    ╎    ╎ 5      @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎   19     …atistics.jl:174; mean
     ╎    ╎    ╎    ╎    ╎    19     …atistics.jl:174; #mean#2
     ╎    ╎    ╎    ╎    ╎     19     …atistics.jl:187; _mean(f::typeof(identit…
     ╎    ╎    ╎    ╎    ╎    ╎ 19     …educedim.jl:1011; sum
     ╎    ╎    ╎    ╎    ╎    ╎  19     …educedim.jl:1011; #sum#829
     ╎    ╎    ╎    ╎    ╎    ╎   19     …educedim.jl:1015; _sum
     ╎    ╎    ╎    ╎    ╎    ╎    19     …educedim.jl:1015; #_sum#831
     ╎    ╎    ╎    ╎    ╎    ╎     19     …educedim.jl:357; mapreduce
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 19     …educedim.jl:357; #mapreduce#821
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  19     …educedim.jl:365; _mapreduce_dim
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   19     …e/reduce.jl:453; _mapreduce
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    19     …e/reduce.jl:175; mapfoldl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     19     …e/reduce.jl:175; #mapfoldl#298
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 19     …e/reduce.jl:44; mapfoldl_impl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  19     …e/reduce.jl:48; foldl_impl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   19     …e/reduce.jl:60; _foldl_impl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    14     …actarray.jl:1215; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     14     …se/range.jl:901; iterate
   14╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 14     …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    5      …actarray.jl:1217; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     5      …subarray.jl:290; getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 5      …subarray.jl:268; reindex
    5╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +2 5      …sentials.jl:13; getindex
    7╎    ╎    ╎    ╎    ╎  8      …pse_comp.jl:124; eclipse_compute_quantities…
    1╎    ╎    ╎    ╎    ╎   1      …sentials.jl:14; getindex(A::Array{Float64,…
   18╎    ╎    ╎    ╎    ╎  19     …pse_comp.jl:125; eclipse_compute_quantities…
    1╎    ╎    ╎    ╎    ╎   1      …se/array.jl:1024; setindex!
     ╎    ╎    ╎    ╎    ╎  19778  …pse_comp.jl:129; eclipse_compute_quantities…
     ╎    ╎    ╎    ╎    ╎   19778  …actarray.jl:3285; map
    8╎    ╎    ╎    ╎    ╎    19778  …se/array.jl:763; collect_similar
     ╎    ╎    ╎    ╎    ╎     154    …se/array.jl:854; _collect(c::Matrix{Floa…
     ╎    ╎    ╎    ╎    ╎    ╎ 154    …enerator.jl:47; iterate
     ╎    ╎    ╎    ╎    ╎    ╎  154    …pse_comp.jl:129; #285
    1╎    ╎    ╎    ╎    ╎    ╎   1      …ase/math.jl:1251; pow_body(x::Float64…
    1╎    ╎    ╎    ╎    ╎    ╎   1      …_physics.jl:0; quad_limb_darkening_ec…
   82╎    ╎    ╎    ╎    ╎    ╎   120    …_physics.jl:29; quad_limb_darkening_e…
     ╎    ╎    ╎    ╎    ╎    ╎    35     …educedim.jl:1153; findmin(f::Functio…
     ╎    ╎    ╎    ╎    ╎    ╎     35     …educedim.jl:1153; #findmin#882
   24╎    ╎    ╎    ╎    ╎    ╎    ╎ 35     …e/reduce.jl:968; _findmin
    6╎    ╎    ╎    ╎    ╎    ╎    ╎  11     …e/reduce.jl:175; mapfoldl
    4╎    ╎    ╎    ╎    ╎    ╎    ╎   5      …e/reduce.jl:175; mapfoldl(f::Fun…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …e/reduce.jl:44; mapfoldl_impl(f…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …e/reduce.jl:48; foldl_impl
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …e/reduce.jl:59; _foldl_impl
    3╎    ╎    ╎    ╎    ╎    ╎    3      …se/tuple.jl:31; getindex(t::Tuple, i…
   25╎    ╎    ╎    ╎    ╎    ╎   32     …_physics.jl:31; quad_limb_darkening_e…
    5╎    ╎    ╎    ╎    ╎    ╎    5      …sentials.jl:13; getindex(A::Vector{F…
     ╎    ╎    ╎    ╎    ╎    ╎    2      …intfuncs.jl:338; literal_pow
     ╎    ╎    ╎    ╎    ╎    ╎     2      …ase/math.jl:1248; ^
    2╎    ╎    ╎    ╎    ╎    ╎    ╎ 2      …ase/math.jl:0; pow_body(x::Float64…
   36╎    ╎    ╎    ╎    ╎     44     …se/array.jl:859; _collect(c::Matrix{Floa…
     ╎    ╎    ╎    ╎    ╎    ╎ 8      …se/array.jl:717; _similar_for(c::Matrix…
     ╎    ╎    ╎    ╎    ╎    ╎  8      …actarray.jl:839; similar
     ╎    ╎    ╎    ╎    ╎    ╎   8      …se/array.jl:420; similar
     ╎    ╎    ╎    ╎    ╎    ╎    8      …ase/boot.jl:487; Array
    8╎    ╎    ╎    ╎    ╎    ╎     8      …ase/boot.jl:479; Array
  103╎    ╎    ╎    ╎    ╎     151    …se/array.jl:863; _collect(c::Matrix{Floa…
   48╎    ╎    ╎    ╎    ╎    ╎ 48     …actarray.jl:274; ndims
   25╎    ╎    ╎    ╎    ╎     19421  …se/array.jl:864; _collect(c::Matrix{Floa…
    2╎    ╎    ╎    ╎    ╎    ╎ 2      …se/array.jl:0; collect_to!(dest::Matrix…
     ╎    ╎    ╎    ╎    ╎    ╎ 19394  …se/array.jl:870; collect_to_with_first!…
     ╎    ╎    ╎    ╎    ╎    ╎  19323  …se/array.jl:892; collect_to!(dest::Mat…
     ╎    ╎    ╎    ╎    ╎    ╎   8      …enerator.jl:44; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    8      …se/array.jl:945; iterate
    8╎    ╎    ╎    ╎    ╎    ╎     8      …sentials.jl:13; getindex
     ╎    ╎    ╎    ╎    ╎    ╎   19315  …enerator.jl:47; iterate
  121╎    ╎    ╎    ╎    ╎    ╎    19315  …pse_comp.jl:129; #285
    5╎    ╎    ╎    ╎    ╎    ╎     5      …ase/math.jl:1251; pow_body(x::Float…
   25╎    ╎    ╎    ╎    ╎    ╎     25     …_physics.jl:0; quad_limb_darkening_…
   30╎    ╎    ╎    ╎    ╎    ╎     30     …_physics.jl:26; quad_limb_darkening…
    1╎    ╎    ╎    ╎    ╎    ╎     8      …_physics.jl:27; quad_limb_darkening…
    7╎    ╎    ╎    ╎    ╎    ╎    ╎ 7      …se/float.jl:536; <
 6626╎    ╎    ╎    ╎    ╎    ╎     11975  …_physics.jl:29; quad_limb_darkening…
   64╎    ╎    ╎    ╎    ╎    ╎    ╎ 5084   …educedim.jl:1153; findmin(f::Funct…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  5020   …educedim.jl:1153; #findmin#882
 2385╎    ╎    ╎    ╎    ╎    ╎    ╎   5020   …e/reduce.jl:968; _findmin
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    2      …terators.jl:281; pairs
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     2      …terators.jl:274; pairs
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 2      …/indices.jl:486; LinearIndices
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  2      …actarray.jl:98; axes
    2╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   2      …se/array.jl:191; size
  353╎    ╎    ╎    ╎    ╎    ╎    ╎    2633   …e/reduce.jl:175; mapfoldl
  759╎    ╎    ╎    ╎    ╎    ╎    ╎     2280   …e/reduce.jl:175; mapfoldl(f::F…
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …e/reduce.jl:0; mapfoldl_impl(…
   18╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 18     …e/reduce.jl:42; mapfoldl_impl…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 9      …e/reduce.jl:43; mapfoldl_impl…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  9      …e/reduce.jl:150; _xfadjoint
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   9      …perators.jl:1041; ComposedF…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    9      …perators.jl:1041; #_#103
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     9      …perators.jl:1044; call_co…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 9      …perators.jl:1045; call_co…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +2 9      …perators.jl:1118; Fix1
    9╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +3 9      …e/reduce.jl:96; MappingRF
   12╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1493   …e/reduce.jl:44; mapfoldl_impl…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  1481   …e/reduce.jl:48; foldl_impl
   62╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   62     …e/reduce.jl:0; _foldl_impl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   12     …e/reduce.jl:56; _foldl_impl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    2      …terators.jl:298; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     2      …/indices.jl:520; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 2      …se/range.jl:897; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +2 2      …se/range.jl:672; isempty
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +3 2      …perators.jl:378; >
    2╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +4 2      @Base/int.jl:83; <
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    10     …terators.jl:301; iterate
    2╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     10     …terators.jl:294; _pairs_e…
    8╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 8      …sentials.jl:13; getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …e/reduce.jl:58; _foldl_impl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …e/reduce.jl:100; MappingRF
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …e/reduce.jl:968; #318
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 1      …_physics.jl:29; #52
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +2 1      …se/float.jl:610; abs
  196╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   196    …e/reduce.jl:59; _foldl_impl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   176    …e/reduce.jl:60; _foldl_impl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    100    …terators.jl:298; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     100    …/indices.jl:520; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 100    …se/range.jl:901; iterate
  100╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +2 100    …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    76     …terators.jl:301; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     76     …terators.jl:294; _pairs_e…
   76╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 76     …sentials.jl:13; getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   1034   …e/reduce.jl:62; _foldl_impl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    1034   …e/reduce.jl:100; MappingRF
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     143    …e/reduce.jl:968; #318
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 143    …_physics.jl:29; #52
   88╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +2 88     …se/float.jl:410; -
   55╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +2 55     …se/float.jl:610; abs
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     891    …e/reduce.jl:86; BottomRF
    6╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 891    …e/reduce.jl:969; _rf_find…
  171╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +2 885    …perators.jl:232; isgreater
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +3 702    …se/float.jl:551; isless
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +4 12     …se/float.jl:544; _fpint
   12╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +5 12     …sentials.jl:581; reinterp…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +4 211    …se/float.jl:545; _fpint
  211╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +5 211    …sentials.jl:647; ifelse
  479╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +4 479    @Base/int.jl:83; <
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +3 12     …perators.jl:247; isunorde…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +4 12     …se/float.jl:620; isnan
   12╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +5 12     …se/float.jl:535; !=
  265╎    ╎    ╎    ╎    ╎    ╎    ╎ 265    …se/tuple.jl:31; getindex(t::Tuple,…
 6585╎    ╎    ╎    ╎    ╎    ╎     7151   …_physics.jl:31; quad_limb_darkening…
  138╎    ╎    ╎    ╎    ╎    ╎    ╎ 138    …sentials.jl:13; getindex(A::Vector…
   68╎    ╎    ╎    ╎    ╎    ╎    ╎ 68     …se/float.jl:411; *(x::Float64, y::…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 5      …intfuncs.jl:332; literal_pow
    5╎    ╎    ╎    ╎    ╎    ╎    ╎  5      …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 5      …intfuncs.jl:333; literal_pow
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  5      …perators.jl:587; *
    5╎    ╎    ╎    ╎    ╎    ╎    ╎   5      …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 309    …intfuncs.jl:338; literal_pow
   15╎    ╎    ╎    ╎    ╎    ╎    ╎  309    …ase/math.jl:1248; ^
   66╎    ╎    ╎    ╎    ╎    ╎    ╎   66     …ase/math.jl:0; pow_body(x::Float…
   10╎    ╎    ╎    ╎    ╎    ╎    ╎   10     …ase/math.jl:1251; pow_body(x::Fl…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   3      …ase/math.jl:1254; pow_body(x::Fl…
    3╎    ╎    ╎    ╎    ╎    ╎    ╎    3      …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   2      …ase/math.jl:1255; pow_body(x::Fl…
    2╎    ╎    ╎    ╎    ╎    ╎    ╎    2      @Base/int.jl:83; <
   12╎    ╎    ╎    ╎    ╎    ╎    ╎   29     …ase/math.jl:1262; pow_body(x::Fl…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    17     …perators.jl:378; >
   17╎    ╎    ╎    ╎    ╎    ╎    ╎     17     @Base/int.jl:83; <
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   4      …ase/math.jl:1263; pow_body(x::Fl…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    4      …perators.jl:378; >
    4╎    ╎    ╎    ╎    ╎    ╎    ╎     4      @Base/int.jl:83; <
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   3      …ase/math.jl:1264; pow_body(x::Fl…
    3╎    ╎    ╎    ╎    ╎    ╎    ╎    3      …se/float.jl:414; muladd
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   7      …ase/math.jl:1265; pow_body(x::Fl…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    7      …ase/math.jl:57; two_mul
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     7      …oatfuncs.jl:439; fma
    7╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 7      …oatfuncs.jl:434; fma_llvm
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   6      …ase/math.jl:1266; pow_body(x::Fl…
    6╎    ╎    ╎    ╎    ╎    ╎    ╎    6      …se/float.jl:409; +
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   20     …ase/math.jl:1268; pow_body(x::Fl…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    20     …perators.jl:587; *
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     20     …romotion.jl:423; *
   20╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 20     …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   15     …ase/math.jl:1269; pow_body(x::Fl…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    15     …ase/math.jl:56; two_mul
   15╎    ╎    ╎    ╎    ╎    ╎    ╎     15     …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   16     …ase/math.jl:1270; pow_body(x::Fl…
   16╎    ╎    ╎    ╎    ╎    ╎    ╎    16     …se/float.jl:409; +
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   5      …ase/math.jl:1271; pow_body(x::Fl…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    5      @Base/int.jl:538; >>>
    5╎    ╎    ╎    ╎    ╎    ╎    ╎     5      @Base/int.jl:530; >>>
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   17     …ase/math.jl:1273; pow_body(x::Fl…
   17╎    ╎    ╎    ╎    ╎    ╎    ╎    17     …se/float.jl:414; muladd
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   91     …ase/math.jl:1274; pow_body(x::Fl…
   74╎    ╎    ╎    ╎    ╎    ╎    ╎    74     …sentials.jl:647; ifelse
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    17     …se/float.jl:623; isfinite
   17╎    ╎    ╎    ╎    ╎    ╎    ╎     17     …se/float.jl:410; -
   21╎    ╎    ╎    ╎    ╎    ╎    ╎ 41     …perators.jl:587; +(::Float64, ::Fl…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  3      …perators.jl:544; afoldl
    3╎    ╎    ╎    ╎    ╎    ╎    ╎   3      …se/float.jl:409; +
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …perators.jl:545; afoldl
    1╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …se/float.jl:409; +
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  16     …perators.jl:546; afoldl
   16╎    ╎    ╎    ╎    ╎    ╎    ╎   16     …se/float.jl:409; +
   10╎    ╎    ╎    ╎    ╎    ╎  10     …se/array.jl:895; collect_to!(dest::Mat…
     ╎    ╎    ╎    ╎    ╎    ╎  5      …se/array.jl:896; collect_to!(dest::Mat…
    5╎    ╎    ╎    ╎    ╎    ╎   5      …se/array.jl:1021; setindex!
     ╎    ╎    ╎    ╎    ╎    ╎  48     …se/array.jl:897; collect_to!(dest::Mat…
   48╎    ╎    ╎    ╎    ╎    ╎   48     @Base/int.jl:87; +
    8╎    ╎    ╎    ╎    ╎    ╎  8      …_physics.jl:26; quad_limb_darkening_ec…
   10╎    ╎    ╎    ╎    ╎  10     …pse_comp.jl:132; eclipse_compute_quantities…
    5╎    ╎    ╎    ╎    ╎  6      …pse_comp.jl:133; eclipse_compute_quantities…
     ╎    ╎    ╎    ╎    ╎   1      …recision.jl:435; step(r::StepRangeLen{Floa…
     ╎    ╎    ╎    ╎    ╎    1      …recision.jl:265; Number
    1╎    ╎    ╎    ╎    ╎     1      …se/float.jl:409; +
   34╎    ╎    ╎    ╎    ╎  560    …pse_comp.jl:134; eclipse_compute_quantities…
   40╎    ╎    ╎    ╎    ╎   526    …actarray.jl:3313; map(f::Function, A::Base…
    2╎    ╎    ╎    ╎    ╎    2      …se/array.jl:827; collect(itr::Base.Genera…
     ╎    ╎    ╎    ╎    ╎    11     …se/array.jl:834; collect(itr::Base.Genera…
     ╎    ╎    ╎    ╎    ╎     3      …enerator.jl:44; iterate
     ╎    ╎    ╎    ╎    ╎    ╎ 3      …terators.jl:1102; iterate
     ╎    ╎    ╎    ╎    ╎    ╎  3      …terators.jl:1094; _piterate
     ╎    ╎    ╎    ╎    ╎    ╎   3      …se/range.jl:892; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    3      …se/range.jl:894; iterate
     ╎    ╎    ╎    ╎    ╎    ╎     3      …recision.jl:481; unsafe_getindex
    3╎    ╎    ╎    ╎    ╎    ╎    ╎ 3      @Base/int.jl:86; -
     ╎    ╎    ╎    ╎    ╎     8      …enerator.jl:47; iterate
    4╎    ╎    ╎    ╎    ╎    ╎ 4      …ase/math.jl:1246; ^(x::Float64, n::Int6…
     ╎    ╎    ╎    ╎    ╎    ╎ 4      …pse_comp.jl:1; (::GRASS.var"#286#290"{F…
     ╎    ╎    ╎    ╎    ╎    ╎  4      …geometry.jl:31; calc_dA
     ╎    ╎    ╎    ╎    ╎    ╎   2      …ial/trig.jl:30; sin(x::Float64)
    2╎    ╎    ╎    ╎    ╎    ╎    2      …se/float.jl:610; abs
     ╎    ╎    ╎    ╎    ╎    ╎   1      …ial/trig.jl:41; sin(x::Float64)
     ╎    ╎    ╎    ╎    ╎    ╎    1      …rem_pio2.jl:243; rem_pio2_kernel
     ╎    ╎    ╎    ╎    ╎    ╎     1      …rem_pio2.jl:52; cody_waite_2c_pio2
    1╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …se/float.jl:414; muladd
     ╎    ╎    ╎    ╎    ╎    ╎   1      …ial/trig.jl:46; sin(x::Float64)
     ╎    ╎    ╎    ╎    ╎    ╎    1      …ial/trig.jl:139; cos_kernel
     ╎    ╎    ╎    ╎    ╎    ╎     1      …ase/math.jl:186; evalpoly
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …ase/math.jl:187; macro expansion
    1╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …se/float.jl:414; muladd
     ╎    ╎    ╎    ╎    ╎    25     …se/array.jl:839; collect(itr::Base.Genera…
     ╎    ╎    ╎    ╎    ╎     25     …se/array.jl:723; _array_for
     ╎    ╎    ╎    ╎    ╎    ╎ 25     …actarray.jl:876; similar
     ╎    ╎    ╎    ╎    ╎    ╎  25     …actarray.jl:877; similar
     ╎    ╎    ╎    ╎    ╎    ╎   25     …ase/boot.jl:487; Array
   25╎    ╎    ╎    ╎    ╎    ╎    25     …ase/boot.jl:479; Array
     ╎    ╎    ╎    ╎    ╎    448    …se/array.jl:844; collect(itr::Base.Genera…
     ╎    ╎    ╎    ╎    ╎     448    …se/array.jl:870; collect_to_with_first!
     ╎    ╎    ╎    ╎    ╎    ╎ 441    …se/array.jl:892; collect_to!(dest::Matr…
    7╎    ╎    ╎    ╎    ╎    ╎  48     …enerator.jl:44; iterate
     ╎    ╎    ╎    ╎    ╎    ╎   41     …terators.jl:1123; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    38     …terators.jl:1110; _piterate1
     ╎    ╎    ╎    ╎    ╎    ╎     4      …se/range.jl:893; iterate
    4╎    ╎    ╎    ╎    ╎    ╎    ╎ 4      @Base/int.jl:83; <
     ╎    ╎    ╎    ╎    ╎    ╎     34     …se/range.jl:894; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 9      …recision.jl:482; unsafe_getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  9      …romotion.jl:423; *
    9╎    ╎    ╎    ╎    ╎    ╎    ╎   9      …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 17     …recision.jl:483; unsafe_getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  17     …recision.jl:84; add12
   14╎    ╎    ╎    ╎    ╎    ╎    ╎   14     …sentials.jl:647; ifelse
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   3      …perators.jl:378; >
    3╎    ╎    ╎    ╎    ╎    ╎    ╎    3      …se/float.jl:536; <
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 8      …recision.jl:484; unsafe_getindex
    8╎    ╎    ╎    ╎    ╎    ╎    ╎  8      …se/float.jl:409; +
     ╎    ╎    ╎    ╎    ╎    ╎    3      …terators.jl:1114; _piterate1
     ╎    ╎    ╎    ╎    ╎    ╎     3      …terators.jl:1110; _piterate1
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 3      …se/range.jl:894; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  2      …recision.jl:483; unsafe_getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   2      …recision.jl:84; add12
    2╎    ╎    ╎    ╎    ╎    ╎    ╎    2      …sentials.jl:647; ifelse
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …recision.jl:484; unsafe_getindex
    1╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …se/float.jl:409; +
    6╎    ╎    ╎    ╎    ╎    ╎  393    …enerator.jl:47; iterate
    2╎    ╎    ╎    ╎    ╎    ╎   2      …ial/trig.jl:0; sin(x::Float64)
    7╎    ╎    ╎    ╎    ╎    ╎   7      …ial/trig.jl:29; sin(x::Float64)
   27╎    ╎    ╎    ╎    ╎    ╎   27     …pse_comp.jl:-133; (::GRASS.var"#286#2…
   82╎    ╎    ╎    ╎    ╎    ╎   351    …pse_comp.jl:1; (::GRASS.var"#286#290"…
    5╎    ╎    ╎    ╎    ╎    ╎    5      …ase/math.jl:0; calc_dA
    7╎    ╎    ╎    ╎    ╎    ╎    264    …geometry.jl:31; calc_dA
    2╎    ╎    ╎    ╎    ╎    ╎     2      …se/float.jl:410; -
     ╎    ╎    ╎    ╎    ╎    ╎     5      …ase/math.jl:1196; ^
    5╎    ╎    ╎    ╎    ╎    ╎    ╎ 5      …romotion.jl:521; ==
   10╎    ╎    ╎    ╎    ╎    ╎     87     …ase/math.jl:1204; ^
    6╎    ╎    ╎    ╎    ╎    ╎    ╎ 6      …ase/math.jl:0; ^(x::Float64, n::In…
    7╎    ╎    ╎    ╎    ╎    ╎    ╎ 7      …ase/math.jl:1246; ^(x::Float64, n:…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 5      …ase/math.jl:1247; ^(x::Float64, n:…
    5╎    ╎    ╎    ╎    ╎    ╎    ╎  5      …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 50     …ase/math.jl:1248; ^(x::Float64, n:…
    6╎    ╎    ╎    ╎    ╎    ╎    ╎  6      …ase/math.jl:0; pow_body(x::Float6…
    3╎    ╎    ╎    ╎    ╎    ╎    ╎  3      …ase/math.jl:1251; pow_body(x::Flo…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  7      …ase/math.jl:1262; pow_body(x::Flo…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   7      …perators.jl:378; >
    7╎    ╎    ╎    ╎    ╎    ╎    ╎    7      @Base/int.jl:83; <
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  5      …ase/math.jl:1268; pow_body(x::Flo…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   5      …perators.jl:587; *
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    5      …romotion.jl:423; *
    5╎    ╎    ╎    ╎    ╎    ╎    ╎     5      …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  11     …ase/math.jl:1270; pow_body(x::Flo…
   11╎    ╎    ╎    ╎    ╎    ╎    ╎   11     …se/float.jl:409; +
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  4      …ase/math.jl:1273; pow_body(x::Flo…
    4╎    ╎    ╎    ╎    ╎    ╎    ╎   4      …se/float.jl:414; muladd
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  14     …ase/math.jl:1274; pow_body(x::Flo…
   14╎    ╎    ╎    ╎    ╎    ╎    ╎   14     …sentials.jl:647; ifelse
    9╎    ╎    ╎    ╎    ╎    ╎    ╎ 9      …ase/math.jl:1251; pow_body(x::Floa…
     ╎    ╎    ╎    ╎    ╎    ╎     53     …perators.jl:587; *
   25╎    ╎    ╎    ╎    ╎    ╎    ╎ 25     …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 28     …perators.jl:544; afoldl
   28╎    ╎    ╎    ╎    ╎    ╎    ╎  28     …se/float.jl:411; *
    7╎    ╎    ╎    ╎    ╎    ╎     7      …ial/trig.jl:29; sin(x::Float64)
     ╎    ╎    ╎    ╎    ╎    ╎     5      …ial/trig.jl:30; sin(x::Float64)
    5╎    ╎    ╎    ╎    ╎    ╎    ╎ 5      …se/float.jl:610; abs
    6╎    ╎    ╎    ╎    ╎    ╎     6      …ial/trig.jl:31; sin(x::Float64)
     ╎    ╎    ╎    ╎    ╎    ╎     5      …ial/trig.jl:38; sin(x::Float64)
    5╎    ╎    ╎    ╎    ╎    ╎    ╎ 5      …se/float.jl:635; isinf
     ╎    ╎    ╎    ╎    ╎    ╎     31     …ial/trig.jl:41; sin(x::Float64)
    4╎    ╎    ╎    ╎    ╎    ╎    ╎ 4      …se/float.jl:0; rem_pio2_kernel
    5╎    ╎    ╎    ╎    ╎    ╎    ╎ 5      …rem_pio2.jl:0; rem_pio2_kernel
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 2      …rem_pio2.jl:224; rem_pio2_kernel
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  2      …ase/math.jl:1540; poshighword
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   2      …ase/math.jl:1541; poshighword
    2╎    ╎    ╎    ╎    ╎    ╎    ╎    2      @Base/int.jl:347; &
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 5      …rem_pio2.jl:229; rem_pio2_kernel
    5╎    ╎    ╎    ╎    ╎    ╎    ╎  5      …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 13     …rem_pio2.jl:236; rem_pio2_kernel
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  5      …rem_pio2.jl:53; cody_waite_2c_pio2
    5╎    ╎    ╎    ╎    ╎    ╎    ╎   5      …se/float.jl:414; muladd
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  8      …rem_pio2.jl:54; cody_waite_2c_pio2
    8╎    ╎    ╎    ╎    ╎    ╎    ╎   8      …se/float.jl:414; muladd
    1╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …rem_pio2.jl:242; rem_pio2_kernel
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …rem_pio2.jl:243; rem_pio2_kernel
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …rem_pio2.jl:52; cody_waite_2c_pio2
    1╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …se/float.jl:414; muladd
    6╎    ╎    ╎    ╎    ╎    ╎     6      …ial/trig.jl:45; sin(x::Float64)
     ╎    ╎    ╎    ╎    ╎    ╎     34     …ial/trig.jl:46; sin(x::Float64)
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 15     …ial/trig.jl:139; cos_kernel
    2╎    ╎    ╎    ╎    ╎    ╎    ╎  2      …se/float.jl:409; +
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  9      …ase/math.jl:186; evalpoly
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   9      …ase/math.jl:187; macro expansion
    9╎    ╎    ╎    ╎    ╎    ╎    ╎    9      …se/float.jl:414; muladd
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  4      …perators.jl:587; *
    4╎    ╎    ╎    ╎    ╎    ╎    ╎   4      …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 19     …ial/trig.jl:142; cos_kernel
    7╎    ╎    ╎    ╎    ╎    ╎    ╎  7      …se/float.jl:411; *
   12╎    ╎    ╎    ╎    ╎    ╎    ╎  12     …se/float.jl:409; +
     ╎    ╎    ╎    ╎    ╎    ╎     16     …ial/trig.jl:48; sin(x::Float64)
    1╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …se/float.jl:407; -
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 7      …ial/trig.jl:72; sin_kernel
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  3      …ase/math.jl:186; evalpoly
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   3      …ase/math.jl:187; macro expansion
    3╎    ╎    ╎    ╎    ╎    ╎    ╎    3      …se/float.jl:414; muladd
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  4      …perators.jl:587; *
    4╎    ╎    ╎    ╎    ╎    ╎    ╎   4      …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 8      …ial/trig.jl:74; sin_kernel
    3╎    ╎    ╎    ╎    ╎    ╎    ╎  3      …se/float.jl:411; *
    5╎    ╎    ╎    ╎    ╎    ╎    ╎  5      …se/float.jl:410; -
    7╎    ╎    ╎    ╎    ╎    ╎ 7      …pse_comp.jl:1; (::GRASS.var"#286#290"{F…
   24╎    ╎    ╎    ╎    ╎  61     …pse_comp.jl:136; eclipse_compute_quantities…
     ╎    ╎    ╎    ╎    ╎   37     …roadcast.jl:903; materialize(bc::Base.Broa…
     ╎    ╎    ╎    ╎    ╎    36     …roadcast.jl:928; copy
     ╎    ╎    ╎    ╎    ╎     28     …roadcast.jl:956; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎ 28     …roadcast.jl:1003; copyto!
    2╎    ╎    ╎    ╎    ╎    ╎  2      …simdloop.jl:0; macro expansion
    5╎    ╎    ╎    ╎    ╎    ╎  5      …simdloop.jl:75; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎  20     …simdloop.jl:77; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎   20     …roadcast.jl:1004; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎    9      …roadcast.jl:636; getindex
     ╎    ╎    ╎    ╎    ╎    ╎     8      …roadcast.jl:681; _broadcast_getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 8      …roadcast.jl:705; _getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  8      …roadcast.jl:675; _broadcast_getin…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   8      …ensional.jl:696; getindex
    8╎    ╎    ╎    ╎    ╎    ╎    ╎    8      …sentials.jl:14; getindex
     ╎    ╎    ╎    ╎    ╎    ╎     1      …roadcast.jl:682; _broadcast_getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …roadcast.jl:709; _broadcast_getind…
    1╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎    11     …ensional.jl:698; setindex!
   11╎    ╎    ╎    ╎    ╎    ╎     11     …se/array.jl:1024; setindex!
     ╎    ╎    ╎    ╎    ╎    ╎  1      …simdloop.jl:78; macro expansion
    1╎    ╎    ╎    ╎    ╎    ╎   1      @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎     8      …roadcast.jl:223; similar
     ╎    ╎    ╎    ╎    ╎    ╎ 8      …roadcast.jl:224; similar
     ╎    ╎    ╎    ╎    ╎    ╎  8      …actarray.jl:876; similar
     ╎    ╎    ╎    ╎    ╎    ╎   8      …actarray.jl:877; similar
     ╎    ╎    ╎    ╎    ╎    ╎    8      …ase/boot.jl:494; Array
     ╎    ╎    ╎    ╎    ╎    ╎     8      …ase/boot.jl:487; Array
    8╎    ╎    ╎    ╎    ╎    ╎    ╎ 8      …ase/boot.jl:479; Array
    1╎    ╎    ╎    ╎    ╎    1      @Base/int.jl:0; copy
   14╎    ╎    ╎    ╎    ╎  137    …pse_comp.jl:137; eclipse_compute_quantities…
    1╎    ╎    ╎    ╎    ╎   1      …se/array.jl:1024; setindex!(::Matrix{Float…
     ╎    ╎    ╎    ╎    ╎   23     …educedim.jl:1010; sum(a::SubArray{Float64,…
     ╎    ╎    ╎    ╎    ╎    23     …educedim.jl:1010; #sum#828
     ╎    ╎    ╎    ╎    ╎     23     …educedim.jl:1014; _sum
     ╎    ╎    ╎    ╎    ╎    ╎ 23     …educedim.jl:1014; #_sum#830
     ╎    ╎    ╎    ╎    ╎    ╎  23     …educedim.jl:1015; _sum
     ╎    ╎    ╎    ╎    ╎    ╎   23     …educedim.jl:1015; #_sum#831
     ╎    ╎    ╎    ╎    ╎    ╎    23     …educedim.jl:357; mapreduce
     ╎    ╎    ╎    ╎    ╎    ╎     23     …educedim.jl:357; #mapreduce#821
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 23     …educedim.jl:365; _mapreduce_dim
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  23     …e/reduce.jl:453; _mapreduce
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   23     …e/reduce.jl:175; mapfoldl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    23     …e/reduce.jl:175; #mapfoldl#298
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     23     …e/reduce.jl:44; mapfoldl_impl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 23     …e/reduce.jl:48; foldl_impl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  23     …e/reduce.jl:60; _foldl_impl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   21     …actarray.jl:1215; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    21     …se/range.jl:901; iterate
   21╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     21     …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   2      …actarray.jl:1217; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    2      …subarray.jl:290; getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     2      …subarray.jl:268; reindex
    2╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 2      …sentials.jl:13; getindex
     ╎    ╎    ╎    ╎    ╎   3      …subarray.jl:183; view(A::Matrix{Float64}, …
     ╎    ╎    ╎    ╎    ╎    2      …/indices.jl:345; to_indices
     ╎    ╎    ╎    ╎    ╎     2      …ensional.jl:869; to_indices
     ╎    ╎    ╎    ╎    ╎    ╎ 2      …ensional.jl:872; _maybe_linear_logical_…
     ╎    ╎    ╎    ╎    ╎    ╎  2      …ensional.jl:789; LogicalIndex
     ╎    ╎    ╎    ╎    ╎    ╎   2      …ensional.jl:785; LogicalIndex
     ╎    ╎    ╎    ╎    ╎    ╎    2      …educedim.jl:439; count
     ╎    ╎    ╎    ╎    ╎    ╎     2      …educedim.jl:439; #count#824
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 2      …educedim.jl:440; count
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  2      …educedim.jl:440; #count#825
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   2      …bitarray.jl:1454; _count
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    2      …bitarray.jl:1446; bitcount
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …bitarray.jl:1449; #bitcount#353
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …sentials.jl:13; getindex
    1╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …bitarray.jl:1450; #bitcount#353
     ╎    ╎    ╎    ╎    ╎    1      …se/tuple.jl:291; map
     ╎    ╎    ╎    ╎    ╎     1      …subarray.jl:183; #214
     ╎    ╎    ╎    ╎    ╎    ╎ 1      …actarray.jl:1481; unalias
     ╎    ╎    ╎    ╎    ╎    ╎  1      …actarray.jl:1516; mightalias
     ╎    ╎    ╎    ╎    ╎    ╎   1      …actarray.jl:1539; dataids
     ╎    ╎    ╎    ╎    ╎    ╎    1      …flection.jl:611; objectid
    1╎    ╎    ╎    ╎    ╎    ╎     1      …flection.jl:617; _objectid
     ╎    ╎    ╎    ╎    ╎   96     …subarray.jl:186; view(A::Matrix{Float64}, …
     ╎    ╎    ╎    ╎    ╎    70     …subarray.jl:126; _maybe_reshape_parent
     ╎    ╎    ╎    ╎    ╎     70     …pedarray.jl:142; reshape
     ╎    ╎    ╎    ╎    ╎    ╎ 70     …pedarray.jl:111; reshape
     ╎    ╎    ╎    ╎    ╎    ╎  1      …actarray.jl:846; to_shape
    1╎    ╎    ╎    ╎    ╎    ╎   1      …se/tuple.jl:291; map
   69╎    ╎    ╎    ╎    ╎    ╎  69     …pedarray.jl:51; reshape
     ╎    ╎    ╎    ╎    ╎    26     …subarray.jl:223; unsafe_view
     ╎    ╎    ╎    ╎    ╎     26     …subarray.jl:28; SubArray
     ╎    ╎    ╎    ╎    ╎    ╎ 26     …ensional.jl:855; ensure_indexable
    1╎    ╎    ╎    ╎    ╎    ╎  26     …ensional.jl:792; collect
     ╎    ╎    ╎    ╎    ╎    ╎   11     …se/array.jl:839; collect(itr::Base.Ge…
     ╎    ╎    ╎    ╎    ╎    ╎    11     …se/array.jl:723; _array_for
     ╎    ╎    ╎    ╎    ╎    ╎     11     …actarray.jl:876; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 11     …actarray.jl:877; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  11     …ase/boot.jl:486; Array
   11╎    ╎    ╎    ╎    ╎    ╎    ╎   11     …ase/boot.jl:477; Array
     ╎    ╎    ╎    ╎    ╎    ╎   14     …se/array.jl:844; collect(itr::Base.Ge…
     ╎    ╎    ╎    ╎    ╎    ╎    14     …se/array.jl:870; collect_to_with_fir…
    1╎    ╎    ╎    ╎    ╎    ╎     1      …se/array.jl:887; collect_to!(dest::…
     ╎    ╎    ╎    ╎    ╎    ╎     7      …se/array.jl:892; collect_to!(dest::…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 7      …enerator.jl:44; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …ensional.jl:836; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      @Base/int.jl:518; ==
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  2      …ensional.jl:841; iterate
    2╎    ╎    ╎    ╎    ╎    ╎    ╎   2      @Base/int.jl:441; trailing_zeros
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …ensional.jl:842; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …bitarray.jl:121; _blsr
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    1      @Base/int.jl:347; &
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  3      …ensional.jl:843; iterate
    3╎    ╎    ╎    ╎    ╎    ╎    ╎   3      @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎    ╎     3      …se/array.jl:896; collect_to!(dest::…
    3╎    ╎    ╎    ╎    ╎    ╎    ╎ 3      …se/array.jl:1021; setindex!
     ╎    ╎    ╎    ╎    ╎    ╎     3      …se/array.jl:897; collect_to!(dest::…
    3╎    ╎    ╎    ╎    ╎    ╎    ╎ 3      @Base/int.jl:87; +
   35╎    ╎    ╎    ╎    ╎  85     …pse_comp.jl:140; eclipse_compute_quantities…
    4╎    ╎    ╎    ╎    ╎   4      …se/array.jl:1024; setindex!(::Matrix{Float…
     ╎    ╎    ╎    ╎    ╎   3      …subarray.jl:183; view(A::Matrix{Float64}, …
     ╎    ╎    ╎    ╎    ╎    1      …/indices.jl:345; to_indices
     ╎    ╎    ╎    ╎    ╎     1      …ensional.jl:869; to_indices
     ╎    ╎    ╎    ╎    ╎    ╎ 1      …ensional.jl:872; _maybe_linear_logical_…
     ╎    ╎    ╎    ╎    ╎    ╎  1      …ensional.jl:789; LogicalIndex
     ╎    ╎    ╎    ╎    ╎    ╎   1      …ensional.jl:785; LogicalIndex
     ╎    ╎    ╎    ╎    ╎    ╎    1      …educedim.jl:439; count
     ╎    ╎    ╎    ╎    ╎    ╎     1      …educedim.jl:439; #count#824
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …educedim.jl:440; count
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …educedim.jl:440; #count#825
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …bitarray.jl:1454; _count
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …bitarray.jl:1446; bitcount
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …bitarray.jl:1449; #bitcount#353
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎    2      …se/tuple.jl:291; map
     ╎    ╎    ╎    ╎    ╎     2      …subarray.jl:183; #214
     ╎    ╎    ╎    ╎    ╎    ╎ 2      …actarray.jl:1481; unalias
     ╎    ╎    ╎    ╎    ╎    ╎  2      …actarray.jl:1516; mightalias
     ╎    ╎    ╎    ╎    ╎    ╎   2      …actarray.jl:1539; dataids
     ╎    ╎    ╎    ╎    ╎    ╎    2      …flection.jl:611; objectid
    2╎    ╎    ╎    ╎    ╎    ╎     2      …flection.jl:617; _objectid
     ╎    ╎    ╎    ╎    ╎   26     …subarray.jl:186; view(A::Matrix{Float64}, …
     ╎    ╎    ╎    ╎    ╎    1      …subarray.jl:126; _maybe_reshape_parent
     ╎    ╎    ╎    ╎    ╎     1      …pedarray.jl:142; reshape
     ╎    ╎    ╎    ╎    ╎    ╎ 1      …pedarray.jl:111; reshape
    1╎    ╎    ╎    ╎    ╎    ╎  1      …pedarray.jl:51; reshape
     ╎    ╎    ╎    ╎    ╎    25     …subarray.jl:223; unsafe_view
     ╎    ╎    ╎    ╎    ╎     25     …subarray.jl:28; SubArray
     ╎    ╎    ╎    ╎    ╎    ╎ 24     …ensional.jl:855; ensure_indexable
     ╎    ╎    ╎    ╎    ╎    ╎  24     …ensional.jl:792; collect
     ╎    ╎    ╎    ╎    ╎    ╎   11     …se/array.jl:839; collect(itr::Base.Ge…
     ╎    ╎    ╎    ╎    ╎    ╎    11     …se/array.jl:723; _array_for
     ╎    ╎    ╎    ╎    ╎    ╎     11     …actarray.jl:876; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 11     …actarray.jl:877; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  11     …ase/boot.jl:486; Array
   11╎    ╎    ╎    ╎    ╎    ╎    ╎   11     …ase/boot.jl:477; Array
     ╎    ╎    ╎    ╎    ╎    ╎   13     …se/array.jl:844; collect(itr::Base.Ge…
     ╎    ╎    ╎    ╎    ╎    ╎    13     …se/array.jl:870; collect_to_with_fir…
     ╎    ╎    ╎    ╎    ╎    ╎     9      …se/array.jl:892; collect_to!(dest::…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 9      …enerator.jl:44; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …ensional.jl:836; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      @Base/int.jl:518; ==
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  4      …ensional.jl:841; iterate
    4╎    ╎    ╎    ╎    ╎    ╎    ╎   4      @Base/int.jl:441; trailing_zeros
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  2      …ensional.jl:842; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   2      …bitarray.jl:121; _blsr
    2╎    ╎    ╎    ╎    ╎    ╎    ╎    2      @Base/int.jl:347; &
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  2      …ensional.jl:843; iterate
    2╎    ╎    ╎    ╎    ╎    ╎    ╎   2      @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎    ╎     1      …se/array.jl:896; collect_to!(dest::…
    1╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …se/array.jl:1021; setindex!
     ╎    ╎    ╎    ╎    ╎    ╎     3      …se/array.jl:897; collect_to!(dest::…
    3╎    ╎    ╎    ╎    ╎    ╎    ╎ 3      @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎    ╎ 1      …subarray.jl:32; SubArray
    1╎    ╎    ╎    ╎    ╎    ╎  1      …subarray.jl:22; SubArray
     ╎    ╎    ╎    ╎    ╎   17     …atistics.jl:174; mean(A::SubArray{Float64,…
     ╎    ╎    ╎    ╎    ╎    17     …atistics.jl:174; #mean#2
     ╎    ╎    ╎    ╎    ╎     17     …atistics.jl:187; _mean(f::typeof(identit…
     ╎    ╎    ╎    ╎    ╎    ╎ 17     …educedim.jl:1011; sum
     ╎    ╎    ╎    ╎    ╎    ╎  17     …educedim.jl:1011; #sum#829
     ╎    ╎    ╎    ╎    ╎    ╎   17     …educedim.jl:1015; _sum
     ╎    ╎    ╎    ╎    ╎    ╎    17     …educedim.jl:1015; #_sum#831
     ╎    ╎    ╎    ╎    ╎    ╎     17     …educedim.jl:357; mapreduce
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 17     …educedim.jl:357; #mapreduce#821
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  17     …educedim.jl:365; _mapreduce_dim
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   17     …e/reduce.jl:453; _mapreduce
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    17     …e/reduce.jl:175; mapfoldl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     17     …e/reduce.jl:175; #mapfoldl#298
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 17     …e/reduce.jl:44; mapfoldl_impl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  17     …e/reduce.jl:48; foldl_impl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   17     …e/reduce.jl:60; _foldl_impl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    13     …actarray.jl:1215; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     13     …se/range.jl:901; iterate
   13╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 13     …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    4      …actarray.jl:1217; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     4      …subarray.jl:290; getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 4      …subarray.jl:268; reindex
    4╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +2 4      …sentials.jl:13; getindex
    1╎    ╎    ╎    ╎    ╎  48     …pse_comp.jl:141; eclipse_compute_quantities…
    2╎    ╎    ╎    ╎    ╎   2      …se/array.jl:1024; setindex!(::Matrix{Float…
     ╎    ╎    ╎    ╎    ╎   1      …subarray.jl:183; view(A::Matrix{Float64}, …
     ╎    ╎    ╎    ╎    ╎    1      …se/tuple.jl:291; map
     ╎    ╎    ╎    ╎    ╎     1      …subarray.jl:183; #214
     ╎    ╎    ╎    ╎    ╎    ╎ 1      …actarray.jl:1481; unalias
     ╎    ╎    ╎    ╎    ╎    ╎  1      …actarray.jl:1516; mightalias
     ╎    ╎    ╎    ╎    ╎    ╎   1      …actarray.jl:1539; dataids
     ╎    ╎    ╎    ╎    ╎    ╎    1      …flection.jl:611; objectid
    1╎    ╎    ╎    ╎    ╎    ╎     1      …flection.jl:617; _objectid
     ╎    ╎    ╎    ╎    ╎   24     …subarray.jl:186; view(A::Matrix{Float64}, …
     ╎    ╎    ╎    ╎    ╎    24     …subarray.jl:223; unsafe_view
     ╎    ╎    ╎    ╎    ╎     24     …subarray.jl:28; SubArray
     ╎    ╎    ╎    ╎    ╎    ╎ 24     …ensional.jl:855; ensure_indexable
     ╎    ╎    ╎    ╎    ╎    ╎  24     …ensional.jl:792; collect
     ╎    ╎    ╎    ╎    ╎    ╎   4      …se/array.jl:839; collect(itr::Base.Ge…
     ╎    ╎    ╎    ╎    ╎    ╎    4      …se/array.jl:723; _array_for
     ╎    ╎    ╎    ╎    ╎    ╎     4      …actarray.jl:876; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 4      …actarray.jl:877; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  4      …ase/boot.jl:486; Array
    4╎    ╎    ╎    ╎    ╎    ╎    ╎   4      …ase/boot.jl:477; Array
     ╎    ╎    ╎    ╎    ╎    ╎   20     …se/array.jl:844; collect(itr::Base.Ge…
     ╎    ╎    ╎    ╎    ╎    ╎    20     …se/array.jl:870; collect_to_with_fir…
     ╎    ╎    ╎    ╎    ╎    ╎     14     …se/array.jl:892; collect_to!(dest::…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 14     …enerator.jl:44; iterate
    1╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …ensional.jl:836; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  3      …ensional.jl:841; iterate
    3╎    ╎    ╎    ╎    ╎    ╎    ╎   3      @Base/int.jl:441; trailing_zeros
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  3      …ensional.jl:842; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   3      …bitarray.jl:121; _blsr
    3╎    ╎    ╎    ╎    ╎    ╎    ╎    3      @Base/int.jl:347; &
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  7      …ensional.jl:843; iterate
    7╎    ╎    ╎    ╎    ╎    ╎    ╎   7      @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎    ╎     2      …se/array.jl:896; collect_to!(dest::…
    2╎    ╎    ╎    ╎    ╎    ╎    ╎ 2      …se/array.jl:1021; setindex!
     ╎    ╎    ╎    ╎    ╎    ╎     4      …se/array.jl:897; collect_to!(dest::…
    4╎    ╎    ╎    ╎    ╎    ╎    ╎ 4      @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎   20     …atistics.jl:174; mean(A::SubArray{Float64,…
     ╎    ╎    ╎    ╎    ╎    20     …atistics.jl:174; #mean#2
     ╎    ╎    ╎    ╎    ╎     20     …atistics.jl:187; _mean(f::typeof(identit…
     ╎    ╎    ╎    ╎    ╎    ╎ 20     …educedim.jl:1011; sum
     ╎    ╎    ╎    ╎    ╎    ╎  20     …educedim.jl:1011; #sum#829
     ╎    ╎    ╎    ╎    ╎    ╎   20     …educedim.jl:1015; _sum
     ╎    ╎    ╎    ╎    ╎    ╎    20     …educedim.jl:1015; #_sum#831
     ╎    ╎    ╎    ╎    ╎    ╎     20     …educedim.jl:357; mapreduce
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 20     …educedim.jl:357; #mapreduce#821
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  20     …educedim.jl:365; _mapreduce_dim
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   20     …e/reduce.jl:453; _mapreduce
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    20     …e/reduce.jl:175; mapfoldl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     20     …e/reduce.jl:175; #mapfoldl#298
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 20     …e/reduce.jl:44; mapfoldl_impl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  20     …e/reduce.jl:48; foldl_impl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   20     …e/reduce.jl:60; _foldl_impl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    14     …actarray.jl:1215; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     14     …se/range.jl:901; iterate
   14╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 14     …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    6      …actarray.jl:1217; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     6      …subarray.jl:290; getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 6      …subarray.jl:268; reindex
    6╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +2 6      …sentials.jl:13; getindex
   25╎    ╎    ╎    ╎    ╎  61     …pse_comp.jl:142; eclipse_compute_quantities…
    4╎    ╎    ╎    ╎    ╎   4      …se/array.jl:1024; setindex!(::Matrix{Float…
     ╎    ╎    ╎    ╎    ╎   15     …educedim.jl:1010; sum(a::SubArray{Float64,…
     ╎    ╎    ╎    ╎    ╎    15     …educedim.jl:1010; #sum#828
     ╎    ╎    ╎    ╎    ╎     15     …educedim.jl:1014; _sum
     ╎    ╎    ╎    ╎    ╎    ╎ 15     …educedim.jl:1014; #_sum#830
     ╎    ╎    ╎    ╎    ╎    ╎  15     …educedim.jl:1015; _sum
     ╎    ╎    ╎    ╎    ╎    ╎   15     …educedim.jl:1015; #_sum#831
     ╎    ╎    ╎    ╎    ╎    ╎    15     …educedim.jl:357; mapreduce
     ╎    ╎    ╎    ╎    ╎    ╎     15     …educedim.jl:357; #mapreduce#821
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 15     …educedim.jl:365; _mapreduce_dim
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  15     …e/reduce.jl:453; _mapreduce
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   15     …e/reduce.jl:175; mapfoldl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    15     …e/reduce.jl:175; #mapfoldl#298
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     15     …e/reduce.jl:44; mapfoldl_impl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 15     …e/reduce.jl:48; foldl_impl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  15     …e/reduce.jl:60; _foldl_impl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   13     …actarray.jl:1215; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    13     …se/range.jl:901; iterate
   13╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     13     …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   2      …actarray.jl:1217; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    2      …subarray.jl:290; getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     2      …subarray.jl:268; reindex
    2╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 2      …sentials.jl:13; getindex
     ╎    ╎    ╎    ╎    ╎   1      …subarray.jl:183; view(A::Matrix{Float64}, …
     ╎    ╎    ╎    ╎    ╎    1      …se/tuple.jl:291; map
     ╎    ╎    ╎    ╎    ╎     1      …subarray.jl:183; #214
     ╎    ╎    ╎    ╎    ╎    ╎ 1      …actarray.jl:1481; unalias
     ╎    ╎    ╎    ╎    ╎    ╎  1      …actarray.jl:1516; mightalias
     ╎    ╎    ╎    ╎    ╎    ╎   1      …actarray.jl:1539; dataids
     ╎    ╎    ╎    ╎    ╎    ╎    1      …flection.jl:611; objectid
    1╎    ╎    ╎    ╎    ╎    ╎     1      …flection.jl:617; _objectid
     ╎    ╎    ╎    ╎    ╎   16     …subarray.jl:186; view(A::Matrix{Float64}, …
     ╎    ╎    ╎    ╎    ╎    16     …subarray.jl:223; unsafe_view
     ╎    ╎    ╎    ╎    ╎     16     …subarray.jl:28; SubArray
     ╎    ╎    ╎    ╎    ╎    ╎ 16     …ensional.jl:855; ensure_indexable
     ╎    ╎    ╎    ╎    ╎    ╎  16     …ensional.jl:792; collect
     ╎    ╎    ╎    ╎    ╎    ╎   3      …se/array.jl:839; collect(itr::Base.Ge…
     ╎    ╎    ╎    ╎    ╎    ╎    3      …se/array.jl:723; _array_for
     ╎    ╎    ╎    ╎    ╎    ╎     3      …actarray.jl:876; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 3      …actarray.jl:877; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  3      …ase/boot.jl:486; Array
    3╎    ╎    ╎    ╎    ╎    ╎    ╎   3      …ase/boot.jl:477; Array
     ╎    ╎    ╎    ╎    ╎    ╎   13     …se/array.jl:844; collect(itr::Base.Ge…
     ╎    ╎    ╎    ╎    ╎    ╎    13     …se/array.jl:870; collect_to_with_fir…
     ╎    ╎    ╎    ╎    ╎    ╎     6      …se/array.jl:892; collect_to!(dest::…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 6      …enerator.jl:44; iterate
    1╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …ensional.jl:836; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  3      …ensional.jl:841; iterate
    3╎    ╎    ╎    ╎    ╎    ╎    ╎   3      @Base/int.jl:441; trailing_zeros
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …ensional.jl:842; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …bitarray.jl:121; _blsr
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    1      @Base/int.jl:347; &
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …ensional.jl:843; iterate
    1╎    ╎    ╎    ╎    ╎    ╎    ╎   1      @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎    ╎     2      …se/array.jl:896; collect_to!(dest::…
    2╎    ╎    ╎    ╎    ╎    ╎    ╎ 2      …se/array.jl:1021; setindex!
     ╎    ╎    ╎    ╎    ╎    ╎     5      …se/array.jl:897; collect_to!(dest::…
    5╎    ╎    ╎    ╎    ╎    ╎    ╎ 5      @Base/int.jl:87; +
   35╎    ╎    ╎    ╎    ╎  87     …pse_comp.jl:143; eclipse_compute_quantities…
    1╎    ╎    ╎    ╎    ╎   1      …se/array.jl:1024; setindex!(::Matrix{Float…
     ╎    ╎    ╎    ╎    ╎   22     …subarray.jl:186; view(A::Matrix{Float64}, …
     ╎    ╎    ╎    ╎    ╎    1      …subarray.jl:126; _maybe_reshape_parent
     ╎    ╎    ╎    ╎    ╎     1      …pedarray.jl:142; reshape
     ╎    ╎    ╎    ╎    ╎    ╎ 1      …pedarray.jl:111; reshape
    1╎    ╎    ╎    ╎    ╎    ╎  1      …pedarray.jl:51; reshape
     ╎    ╎    ╎    ╎    ╎    21     …subarray.jl:223; unsafe_view
     ╎    ╎    ╎    ╎    ╎     21     …subarray.jl:28; SubArray
     ╎    ╎    ╎    ╎    ╎    ╎ 21     …ensional.jl:855; ensure_indexable
     ╎    ╎    ╎    ╎    ╎    ╎  21     …ensional.jl:792; collect
     ╎    ╎    ╎    ╎    ╎    ╎   5      …se/array.jl:839; collect(itr::Base.Ge…
     ╎    ╎    ╎    ╎    ╎    ╎    5      …se/array.jl:723; _array_for
     ╎    ╎    ╎    ╎    ╎    ╎     5      …actarray.jl:876; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 5      …actarray.jl:877; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  5      …ase/boot.jl:486; Array
    5╎    ╎    ╎    ╎    ╎    ╎    ╎   5      …ase/boot.jl:477; Array
     ╎    ╎    ╎    ╎    ╎    ╎   16     …se/array.jl:844; collect(itr::Base.Ge…
     ╎    ╎    ╎    ╎    ╎    ╎    16     …se/array.jl:870; collect_to_with_fir…
     ╎    ╎    ╎    ╎    ╎    ╎     10     …se/array.jl:892; collect_to!(dest::…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 10     …enerator.jl:44; iterate
    1╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …ensional.jl:0; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  2      …ensional.jl:841; iterate
    2╎    ╎    ╎    ╎    ╎    ╎    ╎   2      @Base/int.jl:441; trailing_zeros
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  2      …ensional.jl:842; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   2      …bitarray.jl:121; _blsr
    2╎    ╎    ╎    ╎    ╎    ╎    ╎    2      @Base/int.jl:347; &
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  5      …ensional.jl:843; iterate
    5╎    ╎    ╎    ╎    ╎    ╎    ╎   5      @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎    ╎     1      …se/array.jl:896; collect_to!(dest::…
    1╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …se/array.jl:1021; setindex!
     ╎    ╎    ╎    ╎    ╎    ╎     5      …se/array.jl:897; collect_to!(dest::…
    5╎    ╎    ╎    ╎    ╎    ╎    ╎ 5      @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎   29     …atistics.jl:174; mean(A::SubArray{Float64,…
     ╎    ╎    ╎    ╎    ╎    29     …atistics.jl:174; #mean#2
     ╎    ╎    ╎    ╎    ╎     29     …atistics.jl:187; _mean(f::typeof(identit…
     ╎    ╎    ╎    ╎    ╎    ╎ 29     …educedim.jl:1011; sum
     ╎    ╎    ╎    ╎    ╎    ╎  29     …educedim.jl:1011; #sum#829
     ╎    ╎    ╎    ╎    ╎    ╎   29     …educedim.jl:1015; _sum
     ╎    ╎    ╎    ╎    ╎    ╎    29     …educedim.jl:1015; #_sum#831
     ╎    ╎    ╎    ╎    ╎    ╎     29     …educedim.jl:357; mapreduce
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 29     …educedim.jl:357; #mapreduce#821
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  29     …educedim.jl:365; _mapreduce_dim
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   29     …e/reduce.jl:453; _mapreduce
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    29     …e/reduce.jl:175; mapfoldl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     29     …e/reduce.jl:175; #mapfoldl#298
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 29     …e/reduce.jl:44; mapfoldl_impl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  29     …e/reduce.jl:48; foldl_impl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   29     …e/reduce.jl:60; _foldl_impl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    21     …actarray.jl:1215; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     21     …se/range.jl:901; iterate
   21╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 21     …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    8      …actarray.jl:1217; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     8      …subarray.jl:290; getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 8      …subarray.jl:268; reindex
    8╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +2 8      …sentials.jl:13; getindex
   39╎    ╎    ╎    ╎    ╎  106    …pse_comp.jl:144; eclipse_compute_quantities…
    2╎    ╎    ╎    ╎    ╎   2      …se/array.jl:1024; setindex!(::Matrix{Float…
     ╎    ╎    ╎    ╎    ╎   1      …subarray.jl:183; view(A::Matrix{Float64}, …
     ╎    ╎    ╎    ╎    ╎    1      …/indices.jl:345; to_indices
     ╎    ╎    ╎    ╎    ╎     1      …ensional.jl:869; to_indices
     ╎    ╎    ╎    ╎    ╎    ╎ 1      …ensional.jl:872; _maybe_linear_logical_…
     ╎    ╎    ╎    ╎    ╎    ╎  1      …ensional.jl:789; LogicalIndex
     ╎    ╎    ╎    ╎    ╎    ╎   1      …ensional.jl:785; LogicalIndex
     ╎    ╎    ╎    ╎    ╎    ╎    1      …educedim.jl:439; count
     ╎    ╎    ╎    ╎    ╎    ╎     1      …educedim.jl:439; #count#824
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …educedim.jl:440; count
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …educedim.jl:440; #count#825
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …bitarray.jl:1454; _count
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …bitarray.jl:1446; bitcount
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …bitarray.jl:1449; #bitcount#353
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      @Base/int.jl:415; count_ones
     ╎    ╎    ╎    ╎    ╎   33     …subarray.jl:186; view(A::Matrix{Float64}, …
     ╎    ╎    ╎    ╎    ╎    2      …subarray.jl:126; _maybe_reshape_parent
     ╎    ╎    ╎    ╎    ╎     2      …pedarray.jl:142; reshape
     ╎    ╎    ╎    ╎    ╎    ╎ 2      …pedarray.jl:111; reshape
    2╎    ╎    ╎    ╎    ╎    ╎  2      …pedarray.jl:51; reshape
     ╎    ╎    ╎    ╎    ╎    31     …subarray.jl:223; unsafe_view
     ╎    ╎    ╎    ╎    ╎     31     …subarray.jl:28; SubArray
     ╎    ╎    ╎    ╎    ╎    ╎ 31     …ensional.jl:855; ensure_indexable
     ╎    ╎    ╎    ╎    ╎    ╎  31     …ensional.jl:792; collect
     ╎    ╎    ╎    ╎    ╎    ╎   16     …se/array.jl:839; collect(itr::Base.Ge…
     ╎    ╎    ╎    ╎    ╎    ╎    16     …se/array.jl:723; _array_for
     ╎    ╎    ╎    ╎    ╎    ╎     16     …actarray.jl:876; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 16     …actarray.jl:877; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  16     …ase/boot.jl:486; Array
   16╎    ╎    ╎    ╎    ╎    ╎    ╎   16     …ase/boot.jl:477; Array
     ╎    ╎    ╎    ╎    ╎    ╎   15     …se/array.jl:844; collect(itr::Base.Ge…
     ╎    ╎    ╎    ╎    ╎    ╎    15     …se/array.jl:870; collect_to_with_fir…
     ╎    ╎    ╎    ╎    ╎    ╎     9      …se/array.jl:892; collect_to!(dest::…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 9      …enerator.jl:44; iterate
    1╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …ensional.jl:0; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …ensional.jl:836; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      @Base/int.jl:518; ==
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  2      …ensional.jl:841; iterate
    2╎    ╎    ╎    ╎    ╎    ╎    ╎   2      @Base/int.jl:441; trailing_zeros
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  2      …ensional.jl:842; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   2      …bitarray.jl:121; _blsr
    2╎    ╎    ╎    ╎    ╎    ╎    ╎    2      @Base/int.jl:347; &
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  3      …ensional.jl:843; iterate
    3╎    ╎    ╎    ╎    ╎    ╎    ╎   3      @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎    ╎     6      …se/array.jl:897; collect_to!(dest::…
    6╎    ╎    ╎    ╎    ╎    ╎    ╎ 6      @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎   31     …atistics.jl:174; mean(A::SubArray{Float64,…
     ╎    ╎    ╎    ╎    ╎    31     …atistics.jl:174; #mean#2
     ╎    ╎    ╎    ╎    ╎     31     …atistics.jl:187; _mean(f::typeof(identit…
     ╎    ╎    ╎    ╎    ╎    ╎ 31     …educedim.jl:1011; sum
     ╎    ╎    ╎    ╎    ╎    ╎  31     …educedim.jl:1011; #sum#829
     ╎    ╎    ╎    ╎    ╎    ╎   31     …educedim.jl:1015; _sum
     ╎    ╎    ╎    ╎    ╎    ╎    31     …educedim.jl:1015; #_sum#831
     ╎    ╎    ╎    ╎    ╎    ╎     31     …educedim.jl:357; mapreduce
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 31     …educedim.jl:357; #mapreduce#821
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  31     …educedim.jl:365; _mapreduce_dim
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   31     …e/reduce.jl:453; _mapreduce
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    31     …e/reduce.jl:175; mapfoldl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     31     …e/reduce.jl:175; #mapfoldl#298
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 31     …e/reduce.jl:44; mapfoldl_impl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  31     …e/reduce.jl:48; foldl_impl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   31     …e/reduce.jl:60; _foldl_impl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    20     …actarray.jl:1215; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     20     …se/range.jl:901; iterate
   20╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 20     …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    11     …actarray.jl:1217; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     11     …subarray.jl:290; getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 11     …subarray.jl:268; reindex
   11╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +2 11     …sentials.jl:13; getindex
   38╎    ╎    ╎    ╎    ╎  145    …pse_comp.jl:145; eclipse_compute_quantities…
    1╎    ╎    ╎    ╎    ╎   1      …se/array.jl:1024; setindex!(::Matrix{Float…
     ╎    ╎    ╎    ╎    ╎   21     …roadcast.jl:903; materialize(bc::Base.Broa…
     ╎    ╎    ╎    ╎    ╎    21     …roadcast.jl:928; copy
     ╎    ╎    ╎    ╎    ╎     17     …roadcast.jl:956; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎ 17     …roadcast.jl:1003; copyto!
    3╎    ╎    ╎    ╎    ╎    ╎  3      …simdloop.jl:0; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎  13     …simdloop.jl:77; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎   13     …roadcast.jl:1004; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎    8      …roadcast.jl:636; getindex
     ╎    ╎    ╎    ╎    ╎    ╎     6      …roadcast.jl:681; _broadcast_getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 6      …roadcast.jl:705; _getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  6      …roadcast.jl:675; _broadcast_getin…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   6      …ensional.jl:696; getindex
    6╎    ╎    ╎    ╎    ╎    ╎    ╎    6      …sentials.jl:14; getindex
     ╎    ╎    ╎    ╎    ╎    ╎     2      …roadcast.jl:682; _broadcast_getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 2      …roadcast.jl:709; _broadcast_getind…
    2╎    ╎    ╎    ╎    ╎    ╎    ╎  2      …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎    5      …ensional.jl:698; setindex!
    5╎    ╎    ╎    ╎    ╎    ╎     5      …se/array.jl:1024; setindex!
    1╎    ╎    ╎    ╎    ╎    ╎  1      …simdloop.jl:84; macro expansion
     ╎    ╎    ╎    ╎    ╎     4      …roadcast.jl:223; similar
     ╎    ╎    ╎    ╎    ╎    ╎ 4      …roadcast.jl:224; similar
     ╎    ╎    ╎    ╎    ╎    ╎  4      …actarray.jl:876; similar
     ╎    ╎    ╎    ╎    ╎    ╎   4      …actarray.jl:877; similar
     ╎    ╎    ╎    ╎    ╎    ╎    4      …ase/boot.jl:494; Array
     ╎    ╎    ╎    ╎    ╎    ╎     4      …ase/boot.jl:487; Array
    4╎    ╎    ╎    ╎    ╎    ╎    ╎ 4      …ase/boot.jl:479; Array
     ╎    ╎    ╎    ╎    ╎   1      …subarray.jl:183; view(A::Matrix{Float64}, …
     ╎    ╎    ╎    ╎    ╎    1      …se/tuple.jl:291; map
     ╎    ╎    ╎    ╎    ╎     1      …subarray.jl:183; #214
     ╎    ╎    ╎    ╎    ╎    ╎ 1      …actarray.jl:1481; unalias
     ╎    ╎    ╎    ╎    ╎    ╎  1      …actarray.jl:1516; mightalias
     ╎    ╎    ╎    ╎    ╎    ╎   1      …actarray.jl:1539; dataids
     ╎    ╎    ╎    ╎    ╎    ╎    1      …flection.jl:611; objectid
    1╎    ╎    ╎    ╎    ╎    ╎     1      …flection.jl:617; _objectid
     ╎    ╎    ╎    ╎    ╎   65     …subarray.jl:186; view(A::Matrix{Float64}, …
     ╎    ╎    ╎    ╎    ╎    45     …subarray.jl:126; _maybe_reshape_parent
     ╎    ╎    ╎    ╎    ╎     45     …pedarray.jl:142; reshape
     ╎    ╎    ╎    ╎    ╎    ╎ 45     …pedarray.jl:111; reshape
   45╎    ╎    ╎    ╎    ╎    ╎  45     …pedarray.jl:51; reshape
     ╎    ╎    ╎    ╎    ╎    20     …subarray.jl:223; unsafe_view
     ╎    ╎    ╎    ╎    ╎     20     …subarray.jl:28; SubArray
     ╎    ╎    ╎    ╎    ╎    ╎ 20     …ensional.jl:855; ensure_indexable
     ╎    ╎    ╎    ╎    ╎    ╎  20     …ensional.jl:792; collect
     ╎    ╎    ╎    ╎    ╎    ╎   7      …se/array.jl:839; collect(itr::Base.Ge…
     ╎    ╎    ╎    ╎    ╎    ╎    7      …se/array.jl:723; _array_for
     ╎    ╎    ╎    ╎    ╎    ╎     7      …actarray.jl:876; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 7      …actarray.jl:877; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  7      …ase/boot.jl:486; Array
    7╎    ╎    ╎    ╎    ╎    ╎    ╎   7      …ase/boot.jl:477; Array
     ╎    ╎    ╎    ╎    ╎    ╎   13     …se/array.jl:844; collect(itr::Base.Ge…
     ╎    ╎    ╎    ╎    ╎    ╎    13     …se/array.jl:870; collect_to_with_fir…
     ╎    ╎    ╎    ╎    ╎    ╎     6      …se/array.jl:892; collect_to!(dest::…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 6      …enerator.jl:44; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  2      …ensional.jl:836; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   2      @Base/int.jl:518; ==
    2╎    ╎    ╎    ╎    ╎    ╎    ╎    2      …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  4      …ensional.jl:842; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   4      …bitarray.jl:121; _blsr
    4╎    ╎    ╎    ╎    ╎    ╎    ╎    4      @Base/int.jl:347; &
     ╎    ╎    ╎    ╎    ╎    ╎     1      …se/array.jl:896; collect_to!(dest::…
    1╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …se/array.jl:1021; setindex!
     ╎    ╎    ╎    ╎    ╎    ╎     6      …se/array.jl:897; collect_to!(dest::…
    6╎    ╎    ╎    ╎    ╎    ╎    ╎ 6      @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎   19     …atistics.jl:174; mean(A::SubArray{Float64,…
     ╎    ╎    ╎    ╎    ╎    19     …atistics.jl:174; #mean#2
     ╎    ╎    ╎    ╎    ╎     19     …atistics.jl:187; _mean(f::typeof(identit…
     ╎    ╎    ╎    ╎    ╎    ╎ 19     …educedim.jl:1011; sum
     ╎    ╎    ╎    ╎    ╎    ╎  19     …educedim.jl:1011; #sum#829
     ╎    ╎    ╎    ╎    ╎    ╎   19     …educedim.jl:1015; _sum
     ╎    ╎    ╎    ╎    ╎    ╎    19     …educedim.jl:1015; #_sum#831
     ╎    ╎    ╎    ╎    ╎    ╎     19     …educedim.jl:357; mapreduce
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 19     …educedim.jl:357; #mapreduce#821
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  19     …educedim.jl:365; _mapreduce_dim
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   19     …e/reduce.jl:453; _mapreduce
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    19     …e/reduce.jl:175; mapfoldl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     19     …e/reduce.jl:175; #mapfoldl#298
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 19     …e/reduce.jl:44; mapfoldl_impl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  19     …e/reduce.jl:48; foldl_impl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   19     …e/reduce.jl:60; _foldl_impl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    17     …actarray.jl:1215; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     17     …se/range.jl:901; iterate
   17╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 17     …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    2      …actarray.jl:1217; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     2      …subarray.jl:290; getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 2      …subarray.jl:268; reindex
    2╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +2 2      …sentials.jl:13; getindex
   78╎    ╎    ╎    ╎    ╎  264    …pse_comp.jl:146; eclipse_compute_quantities…
    1╎    ╎    ╎    ╎    ╎   1      …se/array.jl:1023; setindex!(::Matrix{Float…
    1╎    ╎    ╎    ╎    ╎   1      …se/array.jl:1024; setindex!(::Matrix{Float…
    1╎    ╎    ╎    ╎    ╎   73     …roadcast.jl:903; materialize(bc::Base.Broa…
     ╎    ╎    ╎    ╎    ╎    70     …roadcast.jl:928; copy
     ╎    ╎    ╎    ╎    ╎     54     …roadcast.jl:956; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎ 3      …roadcast.jl:1000; copyto!
     ╎    ╎    ╎    ╎    ╎    ╎  3      …roadcast.jl:983; preprocess
     ╎    ╎    ╎    ╎    ╎    ╎   3      …roadcast.jl:986; preprocess_args
     ╎    ╎    ╎    ╎    ╎    ╎    1      …roadcast.jl:983; preprocess
     ╎    ╎    ╎    ╎    ╎    ╎     1      …roadcast.jl:986; preprocess_args
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …roadcast.jl:987; preprocess_args
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …roadcast.jl:984; preprocess
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …roadcast.jl:676; extrude
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …roadcast.jl:625; newindexer
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …roadcast.jl:626; shapeindexer
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …roadcast.jl:630; _newindexer
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …roadcast.jl:631; _newindexer
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …perators.jl:276; !=
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎    2      …roadcast.jl:987; preprocess_args
     ╎    ╎    ╎    ╎    ╎    ╎     2      …roadcast.jl:984; preprocess
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 2      …roadcast.jl:676; extrude
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  2      …roadcast.jl:625; newindexer
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   2      …roadcast.jl:626; shapeindexer
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    2      …roadcast.jl:630; _newindexer
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     2      …roadcast.jl:631; _newindexer
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 2      …perators.jl:276; !=
    2╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  2      …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎ 51     …roadcast.jl:1003; copyto!
    5╎    ╎    ╎    ╎    ╎    ╎  5      …simdloop.jl:0; macro expansion
    1╎    ╎    ╎    ╎    ╎    ╎  1      …simdloop.jl:75; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎  39     …simdloop.jl:77; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎   39     …roadcast.jl:1004; macro expansion
     ╎    ╎    ╎    ╎    ╎    ╎    28     …roadcast.jl:636; getindex
     ╎    ╎    ╎    ╎    ╎    ╎     22     …roadcast.jl:681; _broadcast_getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 22     …roadcast.jl:705; _getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  4      …roadcast.jl:675; _broadcast_getin…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   4      …ensional.jl:696; getindex
    4╎    ╎    ╎    ╎    ╎    ╎    ╎    4      …sentials.jl:14; getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  16     …roadcast.jl:681; _broadcast_getin…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   16     …roadcast.jl:705; _getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    16     …roadcast.jl:675; _broadcast_get…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     16     …ensional.jl:696; getindex
   16╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 16     …sentials.jl:14; getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …roadcast.jl:682; _broadcast_getin…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …roadcast.jl:709; _broadcast_geti…
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …roadcast.jl:706; _getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …roadcast.jl:675; _broadcast_geti…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …ensional.jl:696; getindex
    1╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …sentials.jl:14; getindex
     ╎    ╎    ╎    ╎    ╎    ╎     6      …roadcast.jl:682; _broadcast_getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 6      …roadcast.jl:709; _broadcast_getind…
    6╎    ╎    ╎    ╎    ╎    ╎    ╎  6      …se/float.jl:411; *
     ╎    ╎    ╎    ╎    ╎    ╎    11     …ensional.jl:698; setindex!
   11╎    ╎    ╎    ╎    ╎    ╎     11     …se/array.jl:1024; setindex!
    6╎    ╎    ╎    ╎    ╎    ╎  6      …simdloop.jl:84; macro expansion
     ╎    ╎    ╎    ╎    ╎     16     …roadcast.jl:223; similar
     ╎    ╎    ╎    ╎    ╎    ╎ 16     …roadcast.jl:224; similar
     ╎    ╎    ╎    ╎    ╎    ╎  16     …actarray.jl:876; similar
     ╎    ╎    ╎    ╎    ╎    ╎   16     …actarray.jl:877; similar
     ╎    ╎    ╎    ╎    ╎    ╎    16     …ase/boot.jl:494; Array
     ╎    ╎    ╎    ╎    ╎    ╎     16     …ase/boot.jl:487; Array
   16╎    ╎    ╎    ╎    ╎    ╎    ╎ 16     …ase/boot.jl:479; Array
     ╎    ╎    ╎    ╎    ╎    2      …roadcast.jl:306; instantiate
    1╎    ╎    ╎    ╎    ╎     1      …roadcast.jl:0; combine_axes
     ╎    ╎    ╎    ╎    ╎     1      …roadcast.jl:524; combine_axes
     ╎    ╎    ╎    ╎    ╎    ╎ 1      …roadcast.jl:543; broadcast_shape
    1╎    ╎    ╎    ╎    ╎    ╎  1      …roadcast.jl:0; _bcs
     ╎    ╎    ╎    ╎    ╎   32     …educedim.jl:1010; sum(a::SubArray{Float64,…
     ╎    ╎    ╎    ╎    ╎    32     …educedim.jl:1010; #sum#828
     ╎    ╎    ╎    ╎    ╎     32     …educedim.jl:1014; _sum
     ╎    ╎    ╎    ╎    ╎    ╎ 32     …educedim.jl:1014; #_sum#830
     ╎    ╎    ╎    ╎    ╎    ╎  32     …educedim.jl:1015; _sum
     ╎    ╎    ╎    ╎    ╎    ╎   32     …educedim.jl:1015; #_sum#831
     ╎    ╎    ╎    ╎    ╎    ╎    32     …educedim.jl:357; mapreduce
     ╎    ╎    ╎    ╎    ╎    ╎     32     …educedim.jl:357; #mapreduce#821
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 32     …educedim.jl:365; _mapreduce_dim
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  32     …e/reduce.jl:453; _mapreduce
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   32     …e/reduce.jl:175; mapfoldl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    32     …e/reduce.jl:175; #mapfoldl#298
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     32     …e/reduce.jl:44; mapfoldl_impl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 32     …e/reduce.jl:48; foldl_impl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  32     …e/reduce.jl:60; _foldl_impl
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   27     …actarray.jl:1215; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    27     …se/range.jl:901; iterate
   27╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     27     …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎   5      …actarray.jl:1217; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎    5      …subarray.jl:290; getindex
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎     5      …subarray.jl:268; reindex
    5╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎  +1 5      …sentials.jl:13; getindex
     ╎    ╎    ╎    ╎    ╎   3      …subarray.jl:183; view(A::Matrix{Float64}, …
     ╎    ╎    ╎    ╎    ╎    1      …/indices.jl:345; to_indices
     ╎    ╎    ╎    ╎    ╎     1      …ensional.jl:869; to_indices
     ╎    ╎    ╎    ╎    ╎    ╎ 1      …ensional.jl:872; _maybe_linear_logical_…
     ╎    ╎    ╎    ╎    ╎    ╎  1      …ensional.jl:789; LogicalIndex
     ╎    ╎    ╎    ╎    ╎    ╎   1      …ensional.jl:785; LogicalIndex
     ╎    ╎    ╎    ╎    ╎    ╎    1      …educedim.jl:439; count
     ╎    ╎    ╎    ╎    ╎    ╎     1      …educedim.jl:439; #count#824
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …educedim.jl:440; count
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …educedim.jl:440; #count#825
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   1      …bitarray.jl:1454; _count
     ╎    ╎    ╎    ╎    ╎    ╎    ╎    1      …bitarray.jl:1446; bitcount
     ╎    ╎    ╎    ╎    ╎    ╎    ╎     1      …bitarray.jl:1449; #bitcount#353
    1╎    ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      @Base/int.jl:415; count_ones
     ╎    ╎    ╎    ╎    ╎    2      …se/tuple.jl:291; map
     ╎    ╎    ╎    ╎    ╎     2      …subarray.jl:183; #214
     ╎    ╎    ╎    ╎    ╎    ╎ 2      …actarray.jl:1481; unalias
     ╎    ╎    ╎    ╎    ╎    ╎  2      …actarray.jl:1516; mightalias
     ╎    ╎    ╎    ╎    ╎    ╎   2      …actarray.jl:1539; dataids
     ╎    ╎    ╎    ╎    ╎    ╎    2      …flection.jl:611; objectid
    2╎    ╎    ╎    ╎    ╎    ╎     2      …flection.jl:617; _objectid
     ╎    ╎    ╎    ╎    ╎   76     …subarray.jl:186; view(A::Matrix{Float64}, …
     ╎    ╎    ╎    ╎    ╎    25     …subarray.jl:126; _maybe_reshape_parent
     ╎    ╎    ╎    ╎    ╎     25     …pedarray.jl:142; reshape
     ╎    ╎    ╎    ╎    ╎    ╎ 25     …pedarray.jl:111; reshape
   25╎    ╎    ╎    ╎    ╎    ╎  25     …pedarray.jl:51; reshape
     ╎    ╎    ╎    ╎    ╎    51     …subarray.jl:223; unsafe_view
     ╎    ╎    ╎    ╎    ╎     51     …subarray.jl:28; SubArray
     ╎    ╎    ╎    ╎    ╎    ╎ 51     …ensional.jl:855; ensure_indexable
     ╎    ╎    ╎    ╎    ╎    ╎  51     …ensional.jl:792; collect
     ╎    ╎    ╎    ╎    ╎    ╎   1      …se/array.jl:834; collect(itr::Base.Ge…
     ╎    ╎    ╎    ╎    ╎    ╎    1      …enerator.jl:44; iterate
     ╎    ╎    ╎    ╎    ╎    ╎     1      …ensional.jl:826; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 1      …ensional.jl:843; iterate
    1╎    ╎    ╎    ╎    ╎    ╎    ╎  1      @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎    ╎   13     …se/array.jl:839; collect(itr::Base.Ge…
     ╎    ╎    ╎    ╎    ╎    ╎    13     …se/array.jl:723; _array_for
     ╎    ╎    ╎    ╎    ╎    ╎     13     …actarray.jl:876; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 13     …actarray.jl:877; similar
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  13     …ase/boot.jl:486; Array
   13╎    ╎    ╎    ╎    ╎    ╎    ╎   13     …ase/boot.jl:477; Array
     ╎    ╎    ╎    ╎    ╎    ╎   37     …se/array.jl:844; collect(itr::Base.Ge…
     ╎    ╎    ╎    ╎    ╎    ╎    37     …se/array.jl:870; collect_to_with_fir…
     ╎    ╎    ╎    ╎    ╎    ╎     17     …se/array.jl:892; collect_to!(dest::…
     ╎    ╎    ╎    ╎    ╎    ╎    ╎ 17     …enerator.jl:44; iterate
    1╎    ╎    ╎    ╎    ╎    ╎    ╎  4      …ensional.jl:836; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   3      @Base/int.jl:518; ==
    3╎    ╎    ╎    ╎    ╎    ╎    ╎    3      …romotion.jl:521; ==
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  1      …ensional.jl:841; iterate
    1╎    ╎    ╎    ╎    ╎    ╎    ╎   1      @Base/int.jl:441; trailing_zeros
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  2      …ensional.jl:842; iterate
     ╎    ╎    ╎    ╎    ╎    ╎    ╎   2      …bitarray.jl:121; _blsr
    2╎    ╎    ╎    ╎    ╎    ╎    ╎    2      @Base/int.jl:347; &
     ╎    ╎    ╎    ╎    ╎    ╎    ╎  10     …ensional.jl:843; iterate
   10╎    ╎    ╎    ╎    ╎    ╎    ╎   10     @Base/int.jl:87; +
     ╎    ╎    ╎    ╎    ╎    ╎     7      …se/array.jl:896; collect_to!(dest::…
    7╎    ╎    ╎    ╎    ╎    ╎    ╎ 7      …se/array.jl:1021; setindex!
     ╎    ╎    ╎    ╎    ╎    ╎     13     …se/array.jl:897; collect_to!(dest::…
   13╎    ╎    ╎    ╎    ╎    ╎    ╎ 13     @Base/int.jl:87; +
    8╎    ╎    ╎    ╎    ╎  10     …pse_comp.jl:148; eclipse_compute_quantities…
    2╎    ╎    ╎    ╎    ╎   2      …se/float.jl:620; isnan(x::Float64)
   15╎    ╎    ╎    ╎    ╎  17     …pse_comp.jl:156; eclipse_compute_quantities…
    1╎    ╎    ╎    ╎    ╎   1      …se/range.jl:899; iterate(r::UnitRange{Int6…
     ╎    ╎    ╎    ╎    ╎   1      …se/range.jl:901; iterate(r::UnitRange{Int6…
    1╎    ╎    ╎    ╎    ╎    1      …romotion.jl:521; ==
    1╎    ╎    ╎    ╎    ╎  1      …geometry.jl:50; pole_vector_grid!(A::Matrix…
     ╎1      @Base/client.jl:489; include(fname::String)
     ╎ 1      @Base/loading.jl:2136; _include(mapexpr::Function, mod::Module, _…
     ╎  1      @Base/loading.jl:2076; include_string(mapexpr::typeof(identity),…
     ╎   1      @Base/boot.jl:385; eval
     ╎    1      …venience_eclipse.jl:10; kwcall(::@NamedTuple{verbose::Bool, u…
     ╎     1      …venience_eclipse.jl:18; synthesize_spectra_eclipse(spec::Spe…
     ╎    ╎ 1      …enience_eclipse.jl:62; synth_Eclipse_cpu(spec::SpecParams{F…
     ╎    ╎  1      …isk_sim_eclipse.jl:1; kwcall(::@NamedTuple{skip_times::Bit…
     ╎    ╎   1      …sk_sim_eclipse.jl:11; disk_sim_eclipse(spec::SpecParams{F…
     ╎    ╎    1      …c/eclipse_comp.jl:62; eclipse_compute_quantities!(disk::…
     ╎    ╎     1      …abstractarray.jl:3313; map(f::Function, A::Base.Iterato…
     ╎    ╎    ╎ 1      @Base/array.jl:844; collect(itr::Base.Generator{Base.It…
     ╎    ╎    ╎  1      @Base/array.jl:870; collect_to_with_first!(dest::Matri…
     ╎    ╎    ╎   1      @Base/array.jl:892; collect_to!(dest::Matrix{Float64}…
     ╎    ╎    ╎    1      …se/generator.jl:47; iterate
     ╎    ╎    ╎     1      …eclipse_comp.jl:62; #284
     ╎    ╎    ╎    ╎ 1      …tar_physics.jl:57; v_scalar(lat::Float64, lon::Fl…
    1╎    ╎    ╎    ╎  1      …/operators.jl:587; *
    3╎3      @Base/math.jl:0; pow_body(x::Float64, n::Int64)
     ╎21     @Base/reduce.jl:438; _mapreduce(f::typeof(LinearAlgebra.norm), op:…
   21╎ 21     @Base/essentials.jl:13; getindex
    1╎1      @Base/reduce.jl:0; mapfoldl_impl(f::Base.var"#318#319"{GRASS.var"#…
     ╎12     …gebra/src/generic.jl:598; norm(itr::Vector{Float64}, p::Int64)
    6╎ 12     …Algebra/src/dense.jl:106; norm2
    6╎  6      …ebra/src/generic.jl:462; generic_norm2(x::Vector{Float64})
Total snapshots: 138959. Utilization: 100% across all threads and tasks. Use the `groupby` kwarg to break down by thread and/or task.
