export
    Atom,
    AtomTable,
    atom_by_idx,
    atom_by_name,
    atoms,
    natoms,
    get_full_name,
    is_bound_to,
    is_geminal,
    is_vicinal

"""
    Atom{T} <: AbstractSystemComponent{T}

Mutable representation of an individual atom in a system.

# Public fields
 - `idx::Int`
 - `number::Int`
 - `element::ElementType`
 - `name::String`
 - `atom_type::String`
 - `r::Vector3{T}`
 - `v::Vector3{T}`
 - `F::Vector3{T}`
 - `formal_charge::Int`
 - `charge::T`
 - `radius::T`

# Private fields
 - `properties::Properties`
 - `flags::Flags`
 - `frame_id::Int`
 - `molecule_idx::MaybeInt`
 - `chain_idx::MaybeInt`
 - `fragment_idx::MaybeInt`

# Constructors
```julia
Atom(
    ac::AbstractAtomContainer{T}
    number::Int,
    element::ElementType;
    # keyword arguments
    name::String = "",
    atom_type::String = "",
    r::Vector3{T} = Vector3{T}(0, 0, 0),
    v::Vector3{T} = Vector3{T}(0, 0, 0),
    F::Vector3{T} = Vector3{T}(0, 0, 0),
    formal_charge::Int = 0,
    charge::T = zero(T),
    radius::T = zero(T),
    properties::Properties = Properties(),
    flags::Flags = Flags(),
    frame_id::Int = 1
    molecule_idx::MaybeInt = nothing,
    chain_idx::MaybeInt = nothing,
    fragment_idx::MaybeInt = nothing,
)
```
Creates a new `Atom{T}` in the given atom container.

```julia
Atom(
    number::Int,
    element::ElementType;
    kwargs...
)
```
Creates a new `Atom{Float32}` in the default system. Supports the same keyword arguments as above.
"""
const Atom{T} = SystemComponent{T, _AtomTableRow{T}}

@inline function Atom(
    sys::System{T},
    number::Int,
    element::ElementType;
    frame_id::Int = 1,
    molecule_idx::MaybeInt = nothing,
    chain_idx::MaybeInt = nothing,
    fragment_idx::MaybeInt = nothing,
    kwargs...
) where T
    idx = _next_idx(sys)
    push!(sys._atoms, idx, number, element;
        frame_id = frame_id,
        molecule_idx = molecule_idx,
        chain_idx = chain_idx,
        fragment_idx = fragment_idx,
        kwargs...
    )
    atom_by_idx(sys, idx)
end

@inline function Atom(
    number::Int,
    element::ElementType;
    kwargs...
)
    Atom(default_system(), number, element; kwargs...)
end

"""
    AtomTable{T} <: AbstractSystemComponentTable{T}

Tables.jl-compatible representation of system atoms (or a subset thereof). Atom tables can be
generated using [`atoms`](@ref) or filtered from other atom tables (via `Base.filter`).

# Public columns
 - `idx::AbstractVector{Int}`
 - `number::AbstractVector{Int}`
 - `element::AbstractVector{ElementType}`
 - `name::AbstractVector{String}`
 - `atom_type::AbstractVector{String}`
 - `r::AbstractVector{Vector3{T}}`
 - `v::AbstractVector{Vector3{T}}`
 - `F::AbstractVector{Vector3{T}}`
 - `formal_charge::AbstractVector{Int}`
 - `charge::AbstractVector{T}`
 - `radius::AbstractVector{T}`

# Private columns
 - `properties::AbstractVector{Properties}`
 - `flags::AbstractVector{Flags}`
 - `frame_id::AbstractVector{Int}`
 - `molecule_idx::AbstractVector{MaybeInt}`
 - `chain_idx::AbstractVector{MaybeInt}`
 - `fragment_idx::AbstractVector{MaybeInt}`
"""
const AtomTable{T} = SystemComponentTable{T, Atom{T}}

@inline function _filter_atoms(f::Function, sys::System{T}) where T
    AtomTable{T}(sys, _filter_idx(f, sys._atoms))
end

@inline _table(sys::System{T}, ::Type{Atom{T}}) where T = sys._atoms

@inline function _hascolumn(::Type{<: Atom}, nm::Symbol)
    nm in _atom_table_cols_set || nm in _atom_table_cols_priv
end

@inline function parent_molecule(atom::Atom) 
    isnothing(atom.molecule_idx) ?
        nothing :
        molecule_by_idx(parent(atom), atom.molecule_idx)
end

@inline function parent_chain(atom::Atom)
    isnothing(atom.chain_idx) ?
        nothing :
        chain_by_idx(atom._sys, atom.chain_idx)
end

@inline function parent_fragment(atom::Atom)
    isnothing(atom.fragment_idx) ?
        nothing :
        fragment_by_idx(parent(atom), atom.fragment_idx)
end

"""
    $(TYPEDSIGNATURES)

Returns the `Atom{T}` associated with the given `idx` in `sys`. Throws a `KeyError` if no such
atom exists.
"""
@inline function atom_by_idx(sys::System{T}, idx::Int) where T
    Atom{T}(sys, _row_by_idx(sys._atoms, idx))
end

"""
    $(TYPEDSIGNATURES)

Returns the first `Atom{T}` associated with the given `name` in `ac`. Returns nothing if no such
atom exists.

# Supported keyword arguments
 - `frame_id::MaybeInt = 1`: \
Any value other than `nothing` limits the result to atoms matching this frame ID.
"""
@inline function atom_by_name(
    ac::AbstractAtomContainer{T},
    name::String;
    frame_id::MaybeInt = 1
) where T
    idx = filter(atom -> atom.name == name, atoms(ac; frame_id = frame_id)).idx
    isempty(idx) ? nothing : atom_by_idx(parent(ac), first(idx))
end

"""
    atoms(::Chain)
    atoms(::Fragment)
    atoms(::Molecule)
    atoms(::System)

Returns an `AtomTable{T}` containing all atoms of the given atom container.

# Supported keyword arguments
 - `frame_id::MaybeInt = 1`
 - `molecule_idx::Union{MaybeInt, Some{Nothing}} = nothing`
 - `chain_idx::Union{MaybeInt, Some{Nothing}} = nothing`
 - `fragment_idx::Union{MaybeInt, Some{Nothing}} = nothing`
All keyword arguments limit the results to atoms matching the given IDs. Keyword arguments set to
`nothing` are ignored. You can use `Some(nothing)` to explicitly filter for ID values of `nothing`.
"""
@inline function atoms(sys::System{T};
    frame_id::MaybeInt = 1,
    molecule_idx::Union{MaybeInt, Some{Nothing}} = nothing,
    chain_idx::Union{MaybeInt, Some{Nothing}} = nothing,
    fragment_idx::Union{MaybeInt, Some{Nothing}} = nothing
) where T
    _filter_atoms(atom ->
        (isnothing(frame_id)       || atom.frame_id == frame_id) &&
        (isnothing(molecule_idx)   || atom.molecule_idx == something(molecule_idx)) &&
        (isnothing(chain_idx)      || atom.chain_idx == something(chain_idx)) &&
        (isnothing(fragment_idx)   || atom.fragment_idx == something(fragment_idx)),
        sys
    )
end

"""
    natoms(::Chain)
    natoms(::Fragment)
    natoms(::Molecule)
    natoms(::System)

Returns the number of atoms in the given atom container.

# Supported keyword arguments
See [`atoms`](@ref)
"""
@inline function natoms(sys::System; kwargs...)
    length(atoms(sys; kwargs...))
end

"""
    bonds(::Atom)

Returns a `BondTable{T}` containing all bonds of the given atom.
"""
@inline function bonds(atom::Atom)
    _filter_bonds(
        bond -> bond.a1 == atom.idx || bond.a2 == atom.idx,
        parent(atom)
    )
end

"""
    nbonds(::Atom)

Returns the number of bonds of the given atom.
"""
@inline function nbonds(atom::Atom)
    length(bonds(atom))
end

"""
    push!(::Fragment{T},   ::Atom{T})
    push!(::Molecule{T},   ::Atom{T})
    push!(::System{T},     ::Atom{T})

Creates a copy of the given atom in the given atom container. The new atom is automatically
assigned a new `idx`.

# Supported keyword arguments
See [`atoms`](@ref)
"""
@inline function Base.push!(sys::System{T}, atom::Atom{T};
    frame_id::Int = 1,
    molecule_idx::MaybeInt = nothing,
    chain_idx::MaybeInt = nothing,
    fragment_idx::MaybeInt = nothing
) where T
    Atom(sys, atom.number, atom.element;
        name = atom.name,
        atom_type = atom.atom_type,
        r = atom.r,
        v = atom.v,
        F = atom.F,
        formal_charge = atom.formal_charge,
        charge = atom.charge,
        radius = atom.radius,
        properties = atom.properties,
        flags = atom.flags,
        frame_id = frame_id,
        molecule_idx = molecule_idx,
        chain_idx = chain_idx,
        fragment_idx = fragment_idx
    )
    sys
end

@enumx FullNameType begin
    # Do not add extensions
    NO_VARIANT_EXTENSIONS = 1
    # Add the residue extensions
    ADD_VARIANT_EXTENSIONS = 2
    # Add the residue ID
    ADD_RESIDUE_ID = 3
    # Add the residue ID and the residue extension
    ADD_VARIANT_EXTENSIONS_AND_ID = 4
end

function get_full_name(
    a::Atom{T},
    type::FullNameType.T = FullNameType.ADD_VARIANT_EXTENSIONS
) where {T<:Real}

    # determine the parent's name
    f = parent_fragment(a)

    parent_name = ""

    if isnothing(f)
        # look for a molecule containing the atom
        m = parent_molecule(a)
        parent_name = strip(m.name)
    else
        # retrieve the fragment name
        parent_name = get_full_name(f, type)
    end

    # retrieve the atom name
    name = strip(a.name)

    # add the parent name only if non-empty
    if !isempty(parent_name)
        name = string(parent_name, ":", name)
    end

    name
end

@inline distance(a1::Atom, a2::Atom) = distance(a1.r, a2.r)

"""
    $(TYPEDSIGNATURES)

Decides if two atoms are bound to each other.
Hydrogen bonds (`has_flag(bond, :TYPE__HYDROGEN)`) are ignored.
"""
function is_bound_to(a1::Atom, a2::Atom)
    s = parent(a1)

    if s != parent(a2)
        return false
    end

    return !isnothing(
        findfirst(
            b::Bond -> 
                ((b.a1 == a1.idx) && (b.a2 == a2.idx)) ||
                ((b.a1 == a2.idx) && (b.a2 == a1.idx)), 
            non_hydrogen_bonds(s)
        )
    )
end

"""
    $(TYPEDSIGNATURES)

    Decides if two atoms are geminal.
    
Two atoms are geminal if they do not share a common bond but both have a
bond to a third atom. For example the two hydrogen atoms in water are geminal. 
Hydrogen bonds (`has_flag(bond, :TYPE__HYDROGEN)`) are ignored.
"""
function is_geminal(a1::Atom, a2::Atom)
    if a1 == a2
        return false
    end

    # an atom is geminal to another, if it is not directly bonded to it...
    is_geminal = !is_bound_to(a1, a2)

    # ...and is bonded to an atom that is bonded to the other atom
    is_geminal && any(map(b -> is_bound_to(get_partner(b, a1), a2), non_hydrogen_bonds(a1)))
end

"""
    $(TYPEDSIGNATURES)

Decides if two atoms are vicinal.

Two atoms are vicinal if they are separated by three bonds (1-4 position).
Hydrogen bonds (`has_flag(bond, :TYPE__HYDROGEN)`) are ignored.
"""
function is_vicinal(a1::Atom, a2::Atom)
    if a1 == a2
        return false
    end

    # an atom is vicinal to another, if it is not directly bonded to it...
    is_vicinal = !is_bound_to(a1, a2)

    # ...and is bonded to an atom that is bonded to an atom that is bonded to this atom
    if is_vicinal
        is_vicinal = false

        for b_1 in non_hydrogen_bonds(a1)
            partner_1 = get_partner(b_1, a1)

            for b_2 in non_hydrogen_bonds(partner_1)
                partner_2 = get_partner(b_2, partner_1)

                if is_bound_to(partner_2, a2)
                    return true
                end
            end
        end
    end

    return false
end
