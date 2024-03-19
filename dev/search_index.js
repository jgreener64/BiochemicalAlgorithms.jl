var documenterSearchIndex = {"docs":
[{"location":"public/system/#System-representation","page":"System representation","title":"System representation","text":"","category":"section"},{"location":"public/system/","page":"System representation","title":"System representation","text":"CurrentModule = BiochemicalAlgorithms","category":"page"},{"location":"public/system/","page":"System representation","title":"System representation","text":"Pages = [\"system.md\"]","category":"page"},{"location":"public/system/#Abstract-types","page":"System representation","title":"Abstract types","text":"","category":"section"},{"location":"public/system/","page":"System representation","title":"System representation","text":"AbstractColumnTable\nAbstractSystemComponentTable\nAbstractSystemComponent\nAbstractAtomContainer","category":"page"},{"location":"public/system/#BiochemicalAlgorithms.AbstractColumnTable","page":"System representation","title":"BiochemicalAlgorithms.AbstractColumnTable","text":"abstract type AbstractColumnTable <: Tables.AbstractColumns\n\nAbstract base type for all Tables.jl-compatible column tables.\n\n\n\n\n\n","category":"type"},{"location":"public/system/#BiochemicalAlgorithms.AbstractSystemComponentTable","page":"System representation","title":"BiochemicalAlgorithms.AbstractSystemComponentTable","text":"abstract type AbstractSystemComponentTable{T<:Real} <: AbstractColumnTable\n\nAbstract base type for all Tables.jl-compatible system component tables.\n\n\n\n\n\n","category":"type"},{"location":"public/system/#BiochemicalAlgorithms.AbstractSystemComponent","page":"System representation","title":"BiochemicalAlgorithms.AbstractSystemComponent","text":"abstract type AbstractSystemComponent{T<:Real}\n\nAbstract base type for all components of a system, including the system itself.\n\n\n\n\n\n","category":"type"},{"location":"public/system/#BiochemicalAlgorithms.AbstractAtomContainer","page":"System representation","title":"BiochemicalAlgorithms.AbstractAtomContainer","text":"abstract type AbstractAtomContainer{T} <: AbstractSystemComponent{T}\n\nAbstract base type for all atom containers.\n\n\n\n\n\n","category":"type"},{"location":"public/system/#Common-functions","page":"System representation","title":"Common functions","text":"","category":"section"},{"location":"public/system/","page":"System representation","title":"System representation","text":"has_property\nget_property\nset_property!\nhas_flag\nset_flag!\nunset_flag!","category":"page"},{"location":"public/system/#BiochemicalAlgorithms.has_property","page":"System representation","title":"BiochemicalAlgorithms.has_property","text":"has_property(\n    ac::AbstractSystemComponent,\n    key::Symbol\n) -> Any\n\n\nReturns a Bool indicating whether the given system component has the given property.\n\n\n\n\n\n","category":"function"},{"location":"public/system/#BiochemicalAlgorithms.get_property","page":"System representation","title":"BiochemicalAlgorithms.get_property","text":"get_property(\n    ac::AbstractSystemComponent,\n    key::Symbol\n) -> Any\n\n\nReturns the property associated with the given key in ac.\n\n\n\n\n\nget_property(\n    ac::AbstractSystemComponent,\n    key::Symbol,\n    default\n) -> Any\n\n\nReturns the property associated with the given key in ac. If no such property exists, returns default.\n\n\n\n\n\n","category":"function"},{"location":"public/system/#BiochemicalAlgorithms.set_property!","page":"System representation","title":"BiochemicalAlgorithms.set_property!","text":"set_property!(\n    ac::AbstractSystemComponent,\n    key::Symbol,\n    value\n) -> Any\n\n\nSets the property associated with the given key in ac to the given value.\n\n\n\n\n\n","category":"function"},{"location":"public/system/#BiochemicalAlgorithms.has_flag","page":"System representation","title":"BiochemicalAlgorithms.has_flag","text":"has_flag(ac::AbstractSystemComponent, flag::Symbol) -> Any\n\n\nReturns a Bool indicating whether the given system component has the given flag.\n\n\n\n\n\n","category":"function"},{"location":"public/system/#BiochemicalAlgorithms.set_flag!","page":"System representation","title":"BiochemicalAlgorithms.set_flag!","text":"set_flag!(ac::AbstractSystemComponent, flag::Symbol) -> Any\n\n\nAdds the given flag to ac.\n\n\n\n\n\n","category":"function"},{"location":"public/system/#BiochemicalAlgorithms.unset_flag!","page":"System representation","title":"BiochemicalAlgorithms.unset_flag!","text":"unset_flag!(\n    ac::AbstractSystemComponent,\n    flag::Symbol\n) -> Any\n\n\nRemoves the given flag from ac.\n\n\n\n\n\n","category":"function"},{"location":"public/system/#Systems","page":"System representation","title":"Systems","text":"","category":"section"},{"location":"public/system/","page":"System representation","title":"System representation","text":"System\ndefault_system\nBase.parent(::System)\nparent_system","category":"page"},{"location":"public/system/#BiochemicalAlgorithms.System","page":"System representation","title":"BiochemicalAlgorithms.System","text":"mutable struct System{T} <: AbstractAtomContainer{T}\n\nMutable representation of a biomolecular system.\n\nFields\n\nname::String\nproperties::Properties\nflags::Flags\n\nConstructors\n\nSystem(name::String = \"\", properties::Properties = Properties(), flags::Flags = Flags())\n\nCreates a new and empty System{Float32}.\n\nSystem{T}(name::String = \"\", properties::Properties = Properties(), flags::Flags = Flags())\n\nCreates a new and empty System{T}.\n\n\n\n\n\n","category":"type"},{"location":"public/system/#BiochemicalAlgorithms.default_system","page":"System representation","title":"BiochemicalAlgorithms.default_system","text":"default_system() -> System{Float32}\n\n\nReturns the global default system.\n\n\n\n\n\n","category":"function"},{"location":"public/system/#Base.parent-Tuple{System}","page":"System representation","title":"Base.parent","text":"parent(::Atom)\nparent(::Bond)\nparent(::Chain)\nparent(::Fragment)\nparent(::Molecule)\nparent(::Nucleotide)\nparent(::Residue)\nparent(::System)\n\nReturns the System{T} containing the given object.\n\n\n\n\n\n","category":"method"},{"location":"public/system/#BiochemicalAlgorithms.parent_system","page":"System representation","title":"BiochemicalAlgorithms.parent_system","text":"parent_system(::Atom)\nparent_system(::Bond)\nparent_system(::Chain)\nparent_system(::Fragment)\nparent_system(::Molecule)\nparent_system(::Nucleotide)\nparent_system(::Residue)\nparent_system(::System)\n\nReturns the System{T} containing the given object. Alias for  Base.parent.\n\n\n\n\n\n","category":"function"},{"location":"public/system/#Atoms","page":"System representation","title":"Atoms","text":"","category":"section"},{"location":"public/system/","page":"System representation","title":"System representation","text":"Atom\nAtomTable\natom_by_idx\natom_by_name\natoms\nis_bound_to\nis_geminal\nis_vicinal\nnatoms\nBase.push!(::System{T}, ::Atom{T}) where T","category":"page"},{"location":"public/system/#BiochemicalAlgorithms.Atom","page":"System representation","title":"BiochemicalAlgorithms.Atom","text":"Atom{T} <: AbstractSystemComponent{T}\n\nMutable representation of an individual atom in a system.\n\nPublic fields\n\nidx::Int\nnumber::Int\nelement::ElementType\nname::String\natom_type::String\nr::Vector3{T}\nv::Vector3{T}\nF::Vector3{T}\nformal_charge::Int\ncharge::T\nradius::T\n\nPrivate fields\n\nproperties::Properties\nflags::Flags\nframe_id::Int\nmolecule_idx::MaybeInt\nchain_idx::MaybeInt\nfragment_idx::MaybeInt\nnucleotide_idx::MaybeInt\nresidue_idx::MaybeInt\n\nConstructors\n\nAtom(\n    ac::AbstractAtomContainer{T}\n    number::Int,\n    element::ElementType;\n    # keyword arguments\n    name::String = \"\",\n    atom_type::String = \"\",\n    r::Vector3{T} = Vector3{T}(0, 0, 0),\n    v::Vector3{T} = Vector3{T}(0, 0, 0),\n    F::Vector3{T} = Vector3{T}(0, 0, 0),\n    formal_charge::Int = 0,\n    charge::T = zero(T),\n    radius::T = zero(T),\n    properties::Properties = Properties(),\n    flags::Flags = Flags(),\n    frame_id::Int = 1\n    molecule_idx::MaybeInt = nothing,\n    chain_idx::MaybeInt = nothing,\n    fragment_idx::MaybeInt = nothing,\n    nucleotide_idx::MaybeInt = nothing,\n    residue_idx::MaybeInt = nothing\n)\n\nCreates a new Atom{T} in the given atom container.\n\nAtom(\n    number::Int,\n    element::ElementType;\n    kwargs...\n)\n\nCreates a new Atom{Float32} in the default system. Supports the same keyword arguments as above.\n\n\n\n\n\n","category":"type"},{"location":"public/system/#BiochemicalAlgorithms.AtomTable","page":"System representation","title":"BiochemicalAlgorithms.AtomTable","text":"AtomTable{T} <: AbstractSystemComponentTable{T}\n\nTables.jl-compatible representation of system atoms (or a subset thereof). Atom tables can be generated using atoms or filtered from other atom tables (via Base.filter).\n\nPublic columns\n\nidx::AbstractVector{Int}\nnumber::AbstractVector{Int}\nelement::AbstractVector{ElementType}\nname::AbstractVector{String}\natom_type::AbstractVector{String}\nr::AbstractVector{Vector3{T}}\nv::AbstractVector{Vector3{T}}\nF::AbstractVector{Vector3{T}}\nformal_charge::AbstractVector{Int}\ncharge::AbstractVector{T}\nradius::AbstractVector{T}\n\nPrivate columns\n\nproperties::AbstractVector{Properties}\nflags::AbstractVector{Flags}\nframe_id::AbstractVector{Int}\nmolecule_idx::AbstractVector{MaybeInt}\nchain_idx::AbstractVector{MaybeInt}\nfragment_idx::AbstractVector{MaybeInt}\nnucleotide_idx::AbstractVector{MaybeInt}\nresidue_idx::AbstractVector{MaybeInt}\n\n\n\n\n\n","category":"type"},{"location":"public/system/#BiochemicalAlgorithms.atom_by_idx","page":"System representation","title":"BiochemicalAlgorithms.atom_by_idx","text":"atom_by_idx(sys::System{T}, idx::Int64) -> Atom\n\n\nReturns the Atom{T} associated with the given idx in sys. Throws a KeyError if no such atom exists.\n\n\n\n\n\n","category":"function"},{"location":"public/system/#BiochemicalAlgorithms.atom_by_name","page":"System representation","title":"BiochemicalAlgorithms.atom_by_name","text":"atom_by_name(\n    ac::AbstractAtomContainer{T},\n    name::String;\n    frame_id\n) -> Union{Nothing, Atom}\n\n\nReturns the first Atom{T} associated with the given name in ac. Returns nothing if no such atom exists.\n\nSupported keyword arguments\n\nframe_id::MaybeInt = 1: Any value other than nothing limits the result to atoms matching this frame ID.\n\n\n\n\n\n","category":"function"},{"location":"public/system/#BiochemicalAlgorithms.atoms","page":"System representation","title":"BiochemicalAlgorithms.atoms","text":"atoms(::Chain)\natoms(::Fragment)\natoms(::Molecule)\natoms(::Nucleotide)\natoms(::Residue)\natoms(::System)\n\nReturns an AtomTable{T} containing all atoms of the given atom container.\n\nSupported keyword arguments\n\nframe_id::MaybeInt = 1\nmolecule_idx::Union{MaybeInt, Some{Nothing}} = nothing\nchain_idx::Union{MaybeInt, Some{Nothing}} = nothing\nfragment_idx::Union{MaybeInt, Some{Nothing}} = nothing\nnucleotide_idx::Union{MaybeInt, Some{Nothing}} = nothing\nresidue_idx::Union{MaybeInt, Some{Nothing}} = nothing\n\nAll keyword arguments limit the results to atoms matching the given IDs. Keyword arguments set to nothing are ignored. You can use Some(nothing) to explicitly filter for ID values of nothing.\n\n\n\n\n\natoms(\n    substruct::Substructure{T, A} where A<:AbstractAtomContainer{T};\n    frame_id,\n    molecule_idx,\n    chain_idx,\n    fragment_idx,\n    nucleotide_idx,\n    residue_idx\n) -> SystemComponentTable{T, C} where {T, C<:Atom{T}}\n\n\nReturns an AtomTable for all of the given system's atoms matching the given criteria (value or missing). Fields given as nothing are ignored. The returned table contains all public and private atom fields.\n\n\n\n\n\n","category":"function"},{"location":"public/system/#BiochemicalAlgorithms.is_bound_to","page":"System representation","title":"BiochemicalAlgorithms.is_bound_to","text":"is_bound_to(a1::Atom, a2::Atom) -> Bool\n\n\nDecides if two atoms are bound to each other. Hydrogen bonds (has_flag(bond, :TYPE__HYDROGEN)) are ignored.\n\n\n\n\n\n","category":"function"},{"location":"public/system/#BiochemicalAlgorithms.is_geminal","page":"System representation","title":"BiochemicalAlgorithms.is_geminal","text":"is_geminal(a1::Atom, a2::Atom) -> Union{Missing, Bool}\n\n\nDecides if two atoms are geminal.\n\nTwo atoms are geminal if they do not share a common bond but both have a bond to a third atom. For example the two hydrogen atoms in water are geminal.  Hydrogen bonds (has_flag(bond, :TYPE__HYDROGEN)) are ignored.\n\n\n\n\n\n","category":"function"},{"location":"public/system/#BiochemicalAlgorithms.is_vicinal","page":"System representation","title":"BiochemicalAlgorithms.is_vicinal","text":"is_vicinal(a1::Atom, a2::Atom) -> Bool\n\n\nDecides if two atoms are vicinal.\n\nTwo atoms are vicinal if they are separated by three bonds (1-4 position). Hydrogen bonds (has_flag(bond, :TYPE__HYDROGEN)) are ignored.\n\n\n\n\n\n","category":"function"},{"location":"public/system/#BiochemicalAlgorithms.natoms","page":"System representation","title":"BiochemicalAlgorithms.natoms","text":"natoms(::Chain)\nnatoms(::Fragment)\nnatoms(::Molecule)\nnatoms(::Nucleotide)\nnatoms(::Residue)\nnatoms(::System)\n\nReturns the number of atoms in the given atom container.\n\nSupported keyword arguments\n\nSee atoms\n\n\n\n\n\n","category":"function"},{"location":"public/system/#Base.push!-Union{Tuple{T}, Tuple{System{T}, Atom{T}}} where T","page":"System representation","title":"Base.push!","text":"push!(::Fragment{T},   ::Atom{T})\npush!(::Molecule{T},   ::Atom{T})\npush!(::Nucleotide{T}, ::Atom{T})\npush!(::Residue{T},    ::Atom{T})\npush!(::System{T},     ::Atom{T})\n\nCreates a copy of the given atom in the given atom container. The new atom is automatically assigned a new idx.\n\nSupported keyword arguments\n\nSee atoms\n\n\n\n\n\n","category":"method"},{"location":"public/system/#Bonds","page":"System representation","title":"Bonds","text":"","category":"section"},{"location":"public/system/","page":"System representation","title":"System representation","text":"Bond\nBondTable\nbond_by_idx\nbonds\nnbonds\nBase.push!(::System{T}, ::Bond{T}) where T","category":"page"},{"location":"public/system/#BiochemicalAlgorithms.Bond","page":"System representation","title":"BiochemicalAlgorithms.Bond","text":"Bond{T} <: AbstractAtomContainer{T}\n\nMutable representation of an individual bond in a system.\n\nPublic fields\n\nidx::Int\na1::Int\na2::Int\norder::BondOrderType\n\nPrivate fields\n\nproperties::Properties\nflags::Flags\n\nConstructors\n\nBond(\n    sys::System{T}, \n    a1::Int, \n    a2::Int, \n    order::BondOrderType;\n    # keyword arguments\n    properties::Properties = Properties(),\n    flags::Flags = Flags()\n)\n\nCreates a new Bond{T} in the given system.\n\nBond(\n    a1::Int,\n    a2::Int,\n    order::BondOrderType;\n    # keyword arguments\n    properties::Properties = Properties(),\n    flags::Flags = Flags()\n)\n\nCreates a new Bond{Float32} in the default system.\n\n\n\n\n\n","category":"type"},{"location":"public/system/#BiochemicalAlgorithms.BondTable","page":"System representation","title":"BiochemicalAlgorithms.BondTable","text":"BondTable{T} <: AbstractSystemComponentTable{T}\n\nTables.jl-compatible representation of system bonds (or a subset thereof). Bond tables can be generated using bonds or filtered from other bond tables (via Base.filter).\n\nPublic columns\n\nidx::AbstractVector{Int}\na1::AbstractVector{Int}\na2::AbstractVector{Int}\norder::AbstractVector{BondOrderType}\n\nPrivate columns\n\nproperties::AbstractVector{Properties}\nflags::AbstractVector{Flags}\n\n\n\n\n\n","category":"type"},{"location":"public/system/#BiochemicalAlgorithms.bond_by_idx","page":"System representation","title":"BiochemicalAlgorithms.bond_by_idx","text":"bond_by_idx(sys::System{T}, idx::Int64) -> Bond\n\n\nReturns the Bond{T} associated with the given idx in sys. Throws a KeyError if no such bond exists.\n\n\n\n\n\n","category":"function"},{"location":"public/system/#BiochemicalAlgorithms.bonds","page":"System representation","title":"BiochemicalAlgorithms.bonds","text":"bonds(::Chain)\nbonds(::Fragment)\nbonds(::Molecule)\nbonds(::Nucleotide)\nbonds(::Residue)\nbonds(::System)\n\nReturns a BondTable{T} containing all bonds of the given atom container where at least one associated atom is contained in the same container.\n\nSupported keyword arguments\n\nSee atoms\n\n\n\n\n\nbonds(::Atom)\n\nReturns a BondTable{T} containing all bonds of the given atom.\n\n\n\n\n\n","category":"function"},{"location":"public/system/#BiochemicalAlgorithms.nbonds","page":"System representation","title":"BiochemicalAlgorithms.nbonds","text":"nbonds(::Chain)\nnbonds(::Fragment)\nnbonds(::Molecule)\nnbonds(::Nucleotide)\nnbonds(::Residue)\nnbonds(::System)\n\nReturns the number of bonds in the given atom container where at least one associated atom is contained in the same container.\n\nSupported keyword arguments\n\nSee atoms\n\n\n\n\n\nnbonds(::Atom)\n\nReturns the number of bonds of the given atom.\n\n\n\n\n\n","category":"function"},{"location":"public/system/#Base.push!-Union{Tuple{T}, Tuple{System{T}, Bond{T}}} where T","page":"System representation","title":"Base.push!","text":"push!(::AbstractAtomContainer, ::Bond{T})\n\nCreates a copy of the given bond in the system associated with the given atom container. The new bond is automatically assigned a new idx.\n\n\n\n\n\n","category":"method"},{"location":"public/system/#Molecules","page":"System representation","title":"Molecules","text":"","category":"section"},{"location":"public/system/","page":"System representation","title":"System representation","text":"Molecule\nMoleculeTable\nmolecule_by_idx\nmolecules\nnmolecules\nparent_molecule\nBase.push!(::System{T}, ::Molecule{T}) where T","category":"page"},{"location":"public/system/#BiochemicalAlgorithms.Molecule","page":"System representation","title":"BiochemicalAlgorithms.Molecule","text":"Molecule{T} <: AbstractAtomContainer{T}\n\nMutable representation of an individual molecule in a system.\n\nPublic fields\n\nidx::Int\nname::String\n\nPrivate fields\n\nproperties::Properties\nflags::Flags\n\nConstructors\n\nMolecule(\n    sys::System{T};\n    # keyword arguments\n    name::String = \"\",\n    properties::Properties = Properties(),\n    flags::Flags = Flags()\n)\n\nCreates a new Molecule{T} in the given system.\n\nMolecule(;\n    #keyword arguments\n    name::String = \"\",\n    properties::Properties = Properties(),\n    flags::Flags = Flags()\n)\n\nCreates a new Molecule{Float32} in the default system.\n\n\n\n\n\n","category":"type"},{"location":"public/system/#BiochemicalAlgorithms.MoleculeTable","page":"System representation","title":"BiochemicalAlgorithms.MoleculeTable","text":"MoleculeTable{T} <: AbstractSystemComponentTable{T}\n\nTables.jl-compatible representation of system molecules (or a subset thereof). Molecule tables can be generated using molecules or filtered from other molecule tables (via Base.filter).\n\nPublic columns\n\nidx::AbstractVector{Int}\nname::AbstractVector{String}\n\nPrivate columns\n\nproperties::AbstractVector{Properties}\nflags::AbstractVector{Flags}\n\n\n\n\n\n","category":"type"},{"location":"public/system/#BiochemicalAlgorithms.molecule_by_idx","page":"System representation","title":"BiochemicalAlgorithms.molecule_by_idx","text":"molecule_by_idx(sys::System{T}, idx::Int64) -> Molecule\n\n\nReturns the Molecule{T} associated with the given idx in sys. Throws a KeyError if no such molecule exists.\n\n\n\n\n\n","category":"function"},{"location":"public/system/#BiochemicalAlgorithms.molecules","page":"System representation","title":"BiochemicalAlgorithms.molecules","text":"molecules(sys::System{T}) -> MoleculeTable\n\n\nReturns a MoleculeTable{T} containing all molecules of the given system.\n\n\n\n\n\n","category":"function"},{"location":"public/system/#BiochemicalAlgorithms.nmolecules","page":"System representation","title":"BiochemicalAlgorithms.nmolecules","text":"nmolecules(sys::System) -> Int64\n\n\nReturns the number of molecules in the given system.\n\n\n\n\n\n","category":"function"},{"location":"public/system/#BiochemicalAlgorithms.parent_molecule","page":"System representation","title":"BiochemicalAlgorithms.parent_molecule","text":"parent_molecule(::Atom)\nparent_molecule(::Chain)\nparent_molecule(::Fragment)\nparent_molecule(::Nucleotide)\nparent_molecule(::Residue)\n\nReturns the Molecule{T} containing the given object. Returns nothing if no such molecule exists.\n\n\n\n\n\n","category":"function"},{"location":"public/system/#Base.push!-Union{Tuple{T}, Tuple{System{T}, Molecule{T}}} where T","page":"System representation","title":"Base.push!","text":"push!(::System{T}, ::Molecule{T})\n\nCreates a copy of the given molecule in the given system. The new molecule is automatically assigned a new idx.\n\n\n\n\n\n","category":"method"},{"location":"public/system/#Chains","page":"System representation","title":"Chains","text":"","category":"section"},{"location":"public/system/","page":"System representation","title":"System representation","text":"Chain\nChainTable\nchain_by_idx\nchains\nnchains\nparent_chain\nBase.push!(::Molecule{T}, ::Chain{T}) where T","category":"page"},{"location":"public/system/#BiochemicalAlgorithms.Chain","page":"System representation","title":"BiochemicalAlgorithms.Chain","text":"Chain{T} <: AbstractAtomContainer{T}\n\nMutable representation of an individual chain in a system.\n\nPublic fields\n\nidx::Int\nname::String\n\nPrivate fields\n\nproperties::Properties\nflags::Flags\nmolecule_idx::Int\n\nConstructors\n\nChain(\n    mol::Molecule{T};\n    # keyword arguments\n    name::String = \"\",\n    properties::Properties = Properties(),\n    flags::Flags = Flags()\n)\n\nCreates a new Chain{T} in the given molecule.\n\n\n\n\n\n","category":"type"},{"location":"public/system/#BiochemicalAlgorithms.ChainTable","page":"System representation","title":"BiochemicalAlgorithms.ChainTable","text":"ChainTable{T} <: AbstractSystemComponentTable{T}\n\nTables.jl-compatible representation of system chains (or a subset thereof). Chain tables can be generated using chains or filtered from other chain tables (via Base.filter).\n\nPublic columns\n\nidx::AbstractVector{Int}\nname::AbstractVector{String}\n\nPrivate columns\n\nproperties::AbstractVector{Properties}\nflags::AbstractVector{Flags}\nmolecule_idx::AbstractVector{Int}\n\n\n\n\n\n","category":"type"},{"location":"public/system/#BiochemicalAlgorithms.chain_by_idx","page":"System representation","title":"BiochemicalAlgorithms.chain_by_idx","text":"chain_by_idx(sys::System{T}, idx::Int64) -> Chain\n\n\nReturns the Chain{T} associated with the given idx in sys. Throws a KeyError if no such chain exists.\n\n\n\n\n\n","category":"function"},{"location":"public/system/#BiochemicalAlgorithms.chains","page":"System representation","title":"BiochemicalAlgorithms.chains","text":"chains(::Molecule)\nchains(::System; kwargs...)\n\nReturns a ChainTable{T} containing all chains of the given atom container.\n\nSupported keyword arguments\n\nmolecule_idx::MaybeInt = nothing: Any value other than nothing limits the result to chains belonging to the molecule with the given ID.\n\n\n\n\n\n","category":"function"},{"location":"public/system/#BiochemicalAlgorithms.nchains","page":"System representation","title":"BiochemicalAlgorithms.nchains","text":"nchains(::Molecule)\nnchains(::System; kwargs...)\n\nReturns the number of chains in the given atom container.\n\nSupported keyword arguments\n\nSee chains\n\n\n\n\n\n","category":"function"},{"location":"public/system/#BiochemicalAlgorithms.parent_chain","page":"System representation","title":"BiochemicalAlgorithms.parent_chain","text":"parent_chain(::Atom)\nparent_chain(::Fragment)\nparent_chain(::Nucleotide)\nparent_chain(::Residue)\n\nReturns the Chain{T} containing the given object. Returns nothing if no such chain exists.\n\n\n\n\n\n","category":"function"},{"location":"public/system/#Base.push!-Union{Tuple{T}, Tuple{Molecule{T}, Chain{T}}} where T","page":"System representation","title":"Base.push!","text":"push!(::Molecule{T}, ::Chain{T})\n\nCreates a copy of the given chain in the given molecule. The new chain is automatically assigned a new idx.\n\n\n\n\n\n","category":"method"},{"location":"public/system/#Fragments","page":"System representation","title":"Fragments","text":"","category":"section"},{"location":"public/system/","page":"System representation","title":"System representation","text":"Fragment\nFragmentTable\nfragment_by_idx\nfragments\nnfragments\nparent_fragment\nBase.push!(::Chain{T}, ::Fragment{T}) where T","category":"page"},{"location":"public/system/#BiochemicalAlgorithms.Fragment","page":"System representation","title":"BiochemicalAlgorithms.Fragment","text":"Fragment{T} <: AbstractAtomContainer{T}\n\nMutable representation of an individual fragment in a system.\n\nPublic fields\n\nidx::Int\nnumber::Int\nname::String\n\nPrivate fields\n\nproperties::Properties\nflags::Flags\nmolecule_idx::Int\nchain_idx::Int\n\nConstructors\n\nFragment(\n    chain::Chain{T},\n    number::Int;\n    # keyword arguments\n    name::String = \"\",\n    properties::Properties = Properties(),\n    flags::Flags = Flags()\n)\n\nCreates a new Fragment{T} in the given chain.\n\n\n\n\n\n","category":"type"},{"location":"public/system/#BiochemicalAlgorithms.FragmentTable","page":"System representation","title":"BiochemicalAlgorithms.FragmentTable","text":"FragmentTable{T} <: AbstractSystemComponentTable{T}\n\nTables.jl-compatible representation of system fragments (or a subset thereof). Fragment tables can be generated using fragments or filtered from other fragment tables (via Base.filter).\n\nPublic columns\n\nidx::AbstractVector{Int}\nnumber::AbstractVector{Int}\nname::AbstractVector{String}\n\nPrivate columns\n\nproperties::AbstractVector{Properties}\nflags::AbstractVector{Flags}\nmolecule_idx::AbstractVector{Int}\nchain_idx::AbstractVector{Int}\n\n\n\n\n\n","category":"type"},{"location":"public/system/#BiochemicalAlgorithms.fragment_by_idx","page":"System representation","title":"BiochemicalAlgorithms.fragment_by_idx","text":"fragment_by_idx(sys::System{T}, idx::Int64) -> Fragment\n\n\nReturns the Fragment{T} associated with the given idx in sys. Throws a KeyError if no such fragment exists.\n\n\n\n\n\n","category":"function"},{"location":"public/system/#BiochemicalAlgorithms.fragments","page":"System representation","title":"BiochemicalAlgorithms.fragments","text":"fragments(::Chain)\nfragments(::Molecule)\nfragments(::System)\n\nReturns a FragmentTable{T} containing all fragments of the given atom container.\n\nSupported keyword arguments\n\nmolecule_idx::MaybeInt = nothing\nchain_idx::MaybeInt = nothing\n\nAll keyword arguments limit the results to fragments matching the given IDs. Keyword arguments set to nothing are ignored.\n\n\n\n\n\n","category":"function"},{"location":"public/system/#BiochemicalAlgorithms.nfragments","page":"System representation","title":"BiochemicalAlgorithms.nfragments","text":"nfragments(::Chain)\nnfragments(::Molecule)\nnfragments(::System)\n\nReturns the number of fragments in the given atom container.\n\nSupported keyword arguments\n\nSee fragments\n\n\n\n\n\n","category":"function"},{"location":"public/system/#BiochemicalAlgorithms.parent_fragment","page":"System representation","title":"BiochemicalAlgorithms.parent_fragment","text":"parent_fragment(::Atom)\n\nReturns the Fragment{T} containing the given atom. Returns nothing if no such fragment exists.\n\n\n\n\n\n","category":"function"},{"location":"public/system/#Base.push!-Union{Tuple{T}, Tuple{Chain{T}, Fragment{T}}} where T","page":"System representation","title":"Base.push!","text":"push!(::Chain{T}, ::Fragment{T})\n\nCreates a copy of the given fragment in the given chain. The new fragment is automatically assigned a new idx.\n\n\n\n\n\n","category":"method"},{"location":"public/system/#Nucleotides","page":"System representation","title":"Nucleotides","text":"","category":"section"},{"location":"public/system/","page":"System representation","title":"System representation","text":"Nucleotide\nNucleotideTable\nnnucleotides\nnucleotide_by_idx\nnucleotides\nparent_nucleotide\nBase.push!(::Chain{T}, ::Nucleotide{T}) where T","category":"page"},{"location":"public/system/#BiochemicalAlgorithms.Nucleotide","page":"System representation","title":"BiochemicalAlgorithms.Nucleotide","text":"Nucleotide{T} <: AbstractAtomContainer{T}\n\nMutable representation of an individual nucleotide in a system.\n\nPublic fields\n\nidx::Int\nnumber::Int\nname::String\n\nPrivate fields\n\nproperties::Properties\nflags::Flags\nmolecule_idx::Int\nchain_idx::Int\n\nConstructors\n\nNucleotide(\n    chain::Chain{T},\n    number::Int;\n    # keyword arguments\n    name::String = \"\",\n    properties::Properties = Properties(),\n    flags::Flags = Flags()\n)\n\nCreates a new Nucleotide{T} in the given chain.\n\n\n\n\n\n","category":"type"},{"location":"public/system/#BiochemicalAlgorithms.NucleotideTable","page":"System representation","title":"BiochemicalAlgorithms.NucleotideTable","text":"NucleotideTable{T} <: AbstractSystemComponentTable{T}\n\nTables.jl-compatible representation of system nucleotides (or a subset thereof). Nucleotide tables can be generated using nucleotides or filtered from other nucleotide tables (via Base.filter).\n\nPublic columns\n\nidx::AbstractVector{Int}\nnumber::AbstractVector{Int}\nname::AbstractVector{String}\n\nPrivate columns\n\nproperties::AbstractVector{Properties}\nflags::AbstractVector{Flags}\nmolecule_idx::AbstractVector{Int}\nchain_idx::AbstractVector{Int}\n\n\n\n\n\n","category":"type"},{"location":"public/system/#BiochemicalAlgorithms.nnucleotides","page":"System representation","title":"BiochemicalAlgorithms.nnucleotides","text":"nnucleotides(::Chain)\nnnucleotides(::Molecule)\nnnucleotides(::System)\n\nReturns the number of nucleotides in the given atom container.\n\nSupported keyword arguments\n\nSee nucleotides\n\n\n\n\n\n","category":"function"},{"location":"public/system/#BiochemicalAlgorithms.nucleotide_by_idx","page":"System representation","title":"BiochemicalAlgorithms.nucleotide_by_idx","text":"nucleotide_by_idx(sys::System{T}, idx::Int64) -> Nucleotide\n\n\nReturns the Nucleotide{T} associated with the given idx in sys. Throws a KeyError if no such nucleotide exists.\n\n\n\n\n\n","category":"function"},{"location":"public/system/#BiochemicalAlgorithms.nucleotides","page":"System representation","title":"BiochemicalAlgorithms.nucleotides","text":"nucleotides(::Chain)\nnucleotides(::Molecule)\nnucleotides(::System)\n\nReturns a NucleotideTable{T} containing all nucleotides of the given atom container.\n\nSupported keyword arguments\n\nmolecule_idx::MaybeInt = nothing\nchain_idx::MaybeInt = nothing\n\nAll keyword arguments limit the results to nucleotides matching the given IDs. Keyword arguments set to nothing are ignored.\n\n\n\n\n\n","category":"function"},{"location":"public/system/#BiochemicalAlgorithms.parent_nucleotide","page":"System representation","title":"BiochemicalAlgorithms.parent_nucleotide","text":"parent_nucleotide(::Atom)\n\nReturns the Nucleotide{T} containing the given atom. Returns nothing if no such nucleotide exists.\n\n\n\n\n\n","category":"function"},{"location":"public/system/#Base.push!-Union{Tuple{T}, Tuple{Chain{T}, Nucleotide{T}}} where T","page":"System representation","title":"Base.push!","text":"push!(::Chain{T}, ::Nucleotide{T})\n\nCreates a copy of the given nucleotide in the given chain. The new nucleotide is automatically assigned a new idx.\n\n\n\n\n\n","category":"method"},{"location":"public/system/#Residues","page":"System representation","title":"Residues","text":"","category":"section"},{"location":"public/system/","page":"System representation","title":"System representation","text":"Residue\nResidueTable\nnresidues\nparent_residue\nresidue_by_idx\nresidues\nBase.push!(::Chain{T}, ::Residue{T}) where T","category":"page"},{"location":"public/system/#BiochemicalAlgorithms.Residue","page":"System representation","title":"BiochemicalAlgorithms.Residue","text":"Residue{T} <: AbstractAtomContainer{T}\n\nMutable representation of an individual residue in a system.\n\nPublic fields\n\nidx::Int\nnumber::Int\ntype::AminoAcid\n\nPrivate fields\n\nproperties::Properties\nflags::Flags\nmolecule_idx::Int\nchain_idx::Int\n\nConstructors\n\nResidue(\n    chain::Chain{T},\n    number::Int,\n    type::AminoAcid;\n    # keyword arguments\n    properties::Properties = Properties(),\n    flags::Flags = Flags()\n)\n\nCreates a new Residue{T} in the given chain.\n\n\n\n\n\n","category":"type"},{"location":"public/system/#BiochemicalAlgorithms.ResidueTable","page":"System representation","title":"BiochemicalAlgorithms.ResidueTable","text":"ResidueTable{T} <: AbstractSystemComponentTable{T}\n\nTables.jl-compatible representation of system residues (or a subset thereof). Residue tables can be generated using residues or filtered from other residue tables (via Base.filter).\n\nPublic columns\n\nidx::AbstractVector{Int}\nnumber::AbstractVector{Int}\ntype::AbstractVector{AminoAcid}\n\nPrivate columns\n\nproperties::AbstractVector{Properties}\nflags::AbstractVector{Flags}\nmolecule_idx::AbstractVector{Int}\nchain_idx::AbstractVector{Int}\n\n\n\n\n\n","category":"type"},{"location":"public/system/#BiochemicalAlgorithms.nresidues","page":"System representation","title":"BiochemicalAlgorithms.nresidues","text":"nresidues(::Chain)\nnresidues(::Molecule)\nnresidues(::System)\n\nReturns the number of residues in the given atom container.\n\nSupported keyword arguments\n\nSee residues\n\n\n\n\n\n","category":"function"},{"location":"public/system/#BiochemicalAlgorithms.parent_residue","page":"System representation","title":"BiochemicalAlgorithms.parent_residue","text":"parent_residue(::Atom)\n\nReturns the Residue{T} containing the given atom. Returns nothing if no such residue exists.\n\n\n\n\n\n","category":"function"},{"location":"public/system/#BiochemicalAlgorithms.residue_by_idx","page":"System representation","title":"BiochemicalAlgorithms.residue_by_idx","text":"residue_by_idx(sys::System{T}, idx::Int64) -> Residue\n\n\nReturns the Residue{T} associated with the given idx in sys. Throws a KeyError if no such residue exists.\n\n\n\n\n\n","category":"function"},{"location":"public/system/#BiochemicalAlgorithms.residues","page":"System representation","title":"BiochemicalAlgorithms.residues","text":"residues(::Chain)\nresidues(::Molecule)\nresidues(::System)\n\nReturns a ResidueTable{T} containing all residues of the given atom container.\n\nSupported keyword arguments\n\nmolecule_idx::MaybeInt = nothing\nchain_idx::MaybeInt = nothing\n\nAll keyword arguments limit the results to residues matching the given IDs. Keyword arguments set to nothing are ignored.\n\n\n\n\n\n","category":"function"},{"location":"public/system/#Base.push!-Union{Tuple{T}, Tuple{Chain{T}, Residue{T}}} where T","page":"System representation","title":"Base.push!","text":"push!(::Chain{T}, ::Residue{T})\n\nCreates a copy of the given residue in the given chain. The new residue is automatically assigned a new idx.\n\n\n\n\n\n","category":"method"},{"location":"private/system/#Model","page":"System representation","title":"Model","text":"","category":"section"},{"location":"private/system/","page":"System representation","title":"System representation","text":"CurrentModule = BiochemicalAlgorithms","category":"page"},{"location":"private/system/","page":"System representation","title":"System representation","text":"Pages = [\"system.md\"]","category":"page"},{"location":"private/system/","page":"System representation","title":"System representation","text":"_default_system\n_next_idx","category":"page"},{"location":"private/system/#BiochemicalAlgorithms._default_system","page":"System representation","title":"BiochemicalAlgorithms._default_system","text":"const _default_system\n\nGlobal default system.\n\n\n\n\n\n","category":"constant"},{"location":"private/system/#BiochemicalAlgorithms._next_idx","page":"System representation","title":"BiochemicalAlgorithms._next_idx","text":"_next_idx(sys::System{T}) -> Int64\n\n\nReturns the next available idx for the given system.\n\n\n\n\n\n","category":"function"},{"location":"public/forcefields/#Force-fields","page":"Force fields","title":"Force fields","text":"","category":"section"},{"location":"public/forcefields/","page":"Force fields","title":"Force fields","text":"CurrentModule = BiochemicalAlgorithms","category":"page"},{"location":"public/forcefields/","page":"Force fields","title":"Force fields","text":"Pages = [\"forcefields.md\"]","category":"page"},{"location":"public/forcefields/","page":"Force fields","title":"Force fields","text":"update!\nread_ball_ini_file","category":"page"},{"location":"public/forcefields/#BiochemicalAlgorithms.update!","page":"Force fields","title":"BiochemicalAlgorithms.update!","text":"Update the internal data structures of the force field when the system changes    (e.g., through coordinate updates)\n\nPlease note that changes to the options or the topology require a call to setup!prior to the call toupdate``.\n\n\n\n\n\n","category":"function"},{"location":"public/forcefields/#BiochemicalAlgorithms.read_ball_ini_file","page":"Force fields","title":"BiochemicalAlgorithms.read_ball_ini_file","text":"read_ball_ini_file(path::String; ...) -> BALLIniFile\nread_ball_ini_file(\n    path::String,\n    T;\n    cleanup_keys\n) -> BALLIniFile\n\n\nRead a file in BALL's old Ini file format and return it as a BALLIniFile.\n\nIf cleanup_keys is set to true (the default), keys into the Ini sections are simplified (e.g., ver:version becomes version, key:I becomes I, ...).\n\n\n\n\n\n","category":"function"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = BiochemicalAlgorithms","category":"page"},{"location":"#BiochemicalAlgorithms","page":"Home","title":"BiochemicalAlgorithms","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for BiochemicalAlgorithms.","category":"page"}]
}
