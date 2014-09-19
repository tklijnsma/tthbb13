#this file auto-generates the TTH TTree header

function prefixed(pref::Symbol, t::Type, objs...)
    d = Dict()
    for o::Symbol in objs
        d[symbol(string(pref, "__", o))] = (t, )
    end
    return d
end

#make a set of array branches with dynamic length
#the length is counted by an int n_pref, which is created
function prefixed_dynlength(pref::Symbol, t::Type, objs...; length_branch=nothing)
    d = Dict()
    if length_branch == nothing
        merge!(d, prefixed(:n, Int32, pref))
        veclength = first(keys(d))
    else
        veclength = length_branch
    end
    for o::Symbol in objs
        d[symbol(string(pref, "__", o))] = (t, veclength)
    end
    return d
end

#needs to be hard-coded
#in case we need to store an array per object, we have to give it a fixed size
const M_MAX = 100

const type_map = {
    Float32 => "F",
    Float64 => "D",
    Int32 => "I",
}

const cpp_type_map = {
    Float32 => x->"float $x",
    Float64 => x->"double $x",
    Int32 => x->"int $x",
    Vector{Float32} => x->"float $x[N_MAX]",
    Array{Float32, 2} => x->"float $x[N_MAX][N_MAX]",
    Vector{Int32} => x->"int $x[N_MAX]",
}

const DEF_VAL = {
    Float32 => x->"$x = DEF_VAL_FLOAT",
    Float64 => x->"$x = DEF_VAL_DOUBLE",
    Int32 => x->"$x = DEF_VAL_INT",
    Vector{Float32} => x->"SET_ZERO($x, N_MAX, DEF_VAL_FLOAT)",
    Array{Float32, 2} => x->"SET_ZERO_2($x, N_MAX, DEF_VAL_FLOAT)",
    Vector{Int32} => x->"SET_ZERO($x, N_MAX, DEF_VAL_INT)",
}

branchtype(x) = x
branchtype(t::Tuple) = t[1]

#make a vector branch
function make_branch{T <: Any}(out::IO, bn::Symbol, bt::Type{Vector{T}}, veclength::Symbol)
    write(out, "tree->Branch(\"$bn\", ", "$bn, \"", bn, "[", veclength, "]/", type_map[eltype(bt)], "\");\n")
    ts();ts();write(out, "branch_map[\"$bn\"] = (void*)", bn, ";\n")
end

#make a 2D-array branch
function make_branch{T <: Any}(out::IO, bn::Symbol, bt::Type{Array{T, 2}}, veclength::Symbol)
    write(out, "tree->Branch(\"$bn\", ", "$bn, \"", bn, "[$M_MAX][$M_MAX]/", type_map[eltype(bt)], "\");\n")
    ts();ts();write(out, "branch_map[\"$bn\"] = (void*)", bn, ";\n")
end

function make_branch(out::IO, bn::Symbol, bt::Any)
    write(out, "tree->Branch(\"$bn\", ", "&$bn, \"", bn ,"/", type_map[bt], "\");\n")
    ts();ts();write(out, "branch_map[\"$bn\"] = ", "(void*)&$bn", ";\n")
end


function make_branch_var(out::IO, bn::Symbol, bt::Type)
    write(out, cpp_type_map[bt](bn), ";\n")
end

ts() = write("\t");

function make_class(out::IO, name::Symbol, d::Dict)

    tree = sort(collect(d))
    write(out, "//Autogenerated\n")
    write(out, "#include <TTree.h>\n")
    write(out, "#include <string>\n")
    write(out, "#include <map>\n")
    write(out, "#define N_MAX 500\n")
    write(out, "#define M_MAX $M_MAX\n")
    write(out, "//these are simple 'sentinel values' for uninitialized variables\n")
    write(out, "#define DEF_VAL_FLOAT -9999.0f\n")
    write(out, "#define DEF_VAL_DOUBLE -9999.0d\n")
    write(out, "#define DEF_VAL_INT -9999\n")
    write(out, "#define FLOAT_EPS 0.0000001f\n")
    write(out, "#define DOUBLE_EPS 0.0000001d\n")
    write(out, "constexpr bool is_undef(int x) { return x==DEF_VAL_INT; };\n")
    write(out, "constexpr bool is_undef(float x) { return fabs(x-DEF_VAL_FLOAT) < FLOAT_EPS; };\n")
    write(out, "constexpr bool is_undef(double x) { return fabs(x-DEF_VAL_DOUBLE) < DOUBLE_EPS; };\n")
    write(out, "#define SET_ZERO(x,n,y) for(int i=0;i<n;i++) {x[i]=y;}\n")
    write(out, "#define SET_ZERO_2(x,n,y) for(int i=0;i<n;i++) { for(int j=0;j<n;j++) { x[i][j]=y; } }\n")
    write(out, "class $name {\n")

    write("public:\n")
    ts();write("$name(TTree* _tree);\n")

    ts();write("TTree* tree;\n")
    ts();write("std::map<const std::string, const void*> branch_map;\n")
    for (k, v) in tree
        ts();make_branch_var(out, k, branchtype(v))
    end

    ts();write("void loop_initialize(void) {\n")
    for (k, v) in tree
        ts();ts();write(string(DEF_VAL[branchtype(v)](k), ";\n"))
    end
    ts();write("}\n")

    ts();write("void make_branches(void) {\n")

    #first make Int branches to get dynamic array counters
    for (k, v) in tree
        branchtype(v) != Int32 && continue
        ts();ts();make_branch(out, k, v...)
    end
    for (k, v) in tree
        branchtype(v) == Int32 && continue
        ts();ts();make_branch(out, k, v...)
    end
    ts();write("}\n")

    write(out, "};\n") #end class

end


#####
tree_structure = Dict()

particle_id = [:id]
fourmomentum = [:pt, :eta, :phi, :mass]
fourmomentum_cartesian = [:px, :py, :pz, :e]

#Leptons

merge!(tree_structure,
    prefixed_dynlength(
        :lep, Vector{Float32},
        fourmomentum...,
        :dxy, :dz, :mva,
        :ch_iso, #charged hadron iso
        :puch_iso, #pile-up charged hadron iso
        :ec_iso, #ecal iso
        :hc_iso, #hcal iso
        :p_iso, #particle flow iso
        :ph_iso, #photon iso
        :rel_iso #relative isolation
    )
)

merge!(tree_structure,
    prefixed_dynlength(:lep, Vector{Int32}, particle_id..., :charge, :is_tight, :is_medium, :is_loose, :id_bitmask)
)

merge!(tree_structure,
    prefixed_dynlength(:gen_lep, Vector{Float32}, fourmomentum...; length_branch=:n__lep)
)

merge!(tree_structure,
    prefixed_dynlength(:gen_lep, Vector{Int32}, particle_id..., :status; length_branch=:n__lep)
)

#Jets
merge!(tree_structure,
    prefixed_dynlength(:jet, Vector{Float32},
        fourmomentum...,
        :energy, 
        :bd_csv, #CSV b-discriminator
        :ch_e, #charged hadron energy,
        :ce_e, #charged EM energy,
        :nh_e, #neutral hadron energy,
        :ne_e, #neutral EM energy,
        :mu_e, #muon energy,
        :el_e, #electron energy,
        :ph_e, #photon energy,

    )
)

merge!(tree_structure,
    prefixed_dynlength(:jet, Vector{Int32},
        particle_id...
    )
)


#top tagger jets
merge!(tree_structure,
    prefixed_dynlength(:jet_toptagger, Vector{Float32},
        fourmomentum...,
        :energy, 
        :top_mass, 
        :w_mass,
        :min_mass, #minimum invariant mass pairing 
    )
)
merge!(tree_structure,
    prefixed_dynlength(:jet_toptagger, Vector{Int32},
        :n_sj, #number of subjets 
    )
)

#top jet subjets
merge!(tree_structure,
    prefixed_dynlength(:jet_toptagger_sj, Vector{Float32},
        fourmomentum...,
        :energy, 
    )
)


#jet constituent data
#merge!(tree_structure,
#    prefixed_dynlength(:jet, Array{Float32, 2},
#        :c_pt, :c_eta, :c_phi, :c_id
#    )
#)

for s in [:gen_jet, :gen_jet_parton]
    merge!(tree_structure,
        prefixed_dynlength(s, Vector{Float32}, fourmomentum...; length_branch=:n__jet)
    )
    
    merge!(tree_structure,
        prefixed_dynlength(s, Vector{Int32}, particle_id..., :status; length_branch=:n__jet)
    )
end

#MET
merge!(tree_structure,
    prefixed(:met, Float32, :pt, :phi, :pt__en_up, :pt__en_down)
)

merge!(tree_structure,
    prefixed(:gen_met, Float32, :pt, :phi)
)

#Weights
merge!(tree_structure,
    prefixed(:weight, Float32, :pu, :pu__up, :pu_down, :trigger, :trigger_up, :trigger_down)
)

#Per-event info
merge!(tree_structure,
    prefixed(:event, Int32, :run, :lumi, :id, :json)
)

#sum pt
merge!(tree_structure,
    prefixed(:lhe, Float32, :ht, :n_j)
)

#number of gen-level partons
merge!(tree_structure,
    prefixed(:lhe, Int32, :n_b, :n_c, :n_l, :n_g, :n_e)
)

merge!(tree_structure,
    prefixed(:debug, Float64, :time1r, :time1c)
)

#pv - primary vertices
merge!(tree_structure,
    prefixed(:n, Int32, :pv)
)
merge!(tree_structure,
    prefixed_dynlength(:pvi, Vector{Float32}, :ntrue, :n0)
)
merge!(tree_structure,
    prefixed_dynlength(:pvi, Vector{Int32}, :bx)
)
#####

make_class(STDOUT, :TTHTree, tree_structure)
