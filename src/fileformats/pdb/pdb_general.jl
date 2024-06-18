using Scanf

PDBDefaultOptions = Dict(
    "strict_line_checking" => true,
    "selected_model" => -1,
    "ignore_xplor_pseudo_atoms" => true
)

function parse_element_string(es)
    result = Elements.Unknown

    # handle special cases
    if es == "D"
        result = Elements.H
    elseif es == "X"
        result = Elements.Unknown
    else
        try
            result = parse(ElementType, es)
        catch
            @warn "BiochemicalAlgorithms::PDB::parse_element_string: could not parse element from $(es); returning Unknown"
        end
    end

    return result
end

function handle_record(line, pdb_info, options=PDBDefaultOptions)
    # The PDB format description says: "Each line in the PDB entry file
    # consists of 80 columns."
    if get(options, "strict_line_checking", false)
        if length(line) ≠ 80
            # handle invalid record
            print("Invalid record!")
        end
    end

    # find the record type corresponding to this line
    tag = line[1:6]

    if tag ∉ keys(RECORD_MAP)
        # handle unknown record type
    else
       parse_record(line, RECORD_MAP[tag]; pdb_info=pdb_info, options=options)
    end
end

function ftype(s, T=Float32)
    result = 
        if s == "s"
            String
        elseif s == "ld" || s == "d"
            Int
        elseif s == "f"
            T
        elseif s == "c"
            String
        else
            Nothing
        end
    return result
end

function parse_record(line, record_type; pdb_info::PDBInfo{T}, options=PDBDefaultOptions) where {T}
    fstring = record_type.format_string  
    results = []

    # we start after the tag, which is 6 letters long
    group_start = 7
    for m in eachmatch(r"(\s*)%(-?)(\d+)?(\.\d+)?([[:alpha:]]+)", fstring)
        t = ftype(m.captures[5], T)

        # skip preceeding whitespace
        group_start += length(m.captures[1])

        group_end = group_start + 
            (isnothing(m.captures[3]) ? 1 : parse(Int, m.captures[3]))
        
        s = line[group_start:group_end-1]

        group_start = group_end

        if t == String
            append!(results, [s])
        else
            _, v = Scanf.scanf(s, Scanf.Format("%" * m.captures[5]), t)
            append!(results, [v])
        end
    end
    
    interpret_record(Val(record_type.record_type), record_type.record_type, results...; pdb_info=pdb_info, options=options)
end

function interpret_record(::Val{RECORD_TYPE__HEADER}, tag, classification, deposition_date, id; pdb_info=pdb_info, options=options)
    pdb_info.name = classification
    pdb_info.deposition_date = deposition_date
    pdb_info.id = id
end

function interpret_record(::Val{RECORD_TYPE__TITLE}, tag, continuation, title; pdb_info, options)
    pdb_info.title *= title
end

function interpret_record(::Val{RECORD_TYPE__TER}, tag, serial_number, residue_name, chain_id, 
        residue_number, residue_insertion_code; pdb_info, options)

        # should we skip this model?
    if ((pdb_info.selected_model != -1) && (pdb_info.selected_model != pdb_info.current_model))
        return
    end

    # close the current chain
    pdb_info.current_chain = nothing
end

function interpret_record(
        ::Val{RECORD_TYPE__ATOM},
        tag,
        serial_number,
        atom_name,
        alternate_location_identifier,
        residue_name,
        chain_id,
        residue_sequence_number,
        residue_insertion_code,
        x,
        y,
        z,
        occupancy,
        temperature_factor,
        segment_id,
        element_symbol,
        charge;
        pdb_info::PDBInfo{T},
        options,
        is_hetero_atom=false) where {T}

    # should we skip this model?
    if ((pdb_info.selected_model != -1) && (pdb_info.selected_model != pdb_info.current_model))
        return
    end

    # is this an XPLOR pseudo atom that we should skip?
    if (get(options, "ignore_xplor_pseudo_atoms", true) && x >= 9998.0 && y >= 9998.0 && z >= 9998.0)
        return
    end

    # is this a new chain?
    if (    isnothing(pdb_info.current_chain)
         || chain_id != pdb_info.current_chain.name )
        pdb_info.current_chain = Chain(pdb_info.mol; name=chain_id)
        pdb_info.current_residue = nothing
    end

    # right now, we only read the first alternate location indicator
    if (      strip(alternate_location_identifier) != "" 
          && (alternate_location_identifier != pdb_info.alternate_location_identifier) )
        return
    end

    # is this a new residue?
    if (    isnothing(pdb_info.current_residue) 
         || residue_sequence_number != pdb_info.current_residue.number
         || residue_name != pdb_info.current_residue.name
         || residue_insertion_code != get(pdb_info.current_residue.properties, :insertion_code, "") )

            pdb_info.current_residue = Fragment(
                pdb_info.current_chain, 
                residue_sequence_number; 
                name=residue_name,
                properties=Properties([
                    :is_hetero_fragment    => is_hetero_atom,
                    :insertion_code        => residue_insertion_code,
                    :alternate_location_id => alternate_location_identifier
                ])
            )
    end

    # finally, create the atom
    formal_charge = tryparse(Int, charge)
    formal_charge = isnothing(formal_charge) ? 0 : formal_charge

    Atom(
        pdb_info.current_residue, 
        serial_number, 
        parse_element_string(strip(element_symbol));
        name = atom_name,
        r = Vector3{T}(x, y, z),
        formal_charge = formal_charge,
        properties = Properties([
            :tempfactor            => temperature_factor,
            :occupancy             => occupancy,
            :is_hetero_atom        => is_hetero_atom,
            :insertion_code        => residue_insertion_code,
            :alternate_location_id => alternate_location_identifier
        ]),
        flags = Flags(),
        frame_id = pdb_info.current_model
    )
end


function interpret_record(
        ::Val{RECORD_TYPE__HETATM},
        tag,
        serial_number,
        atom_name,
        alternate_location_identifier,
        residue_name,
        chain_id,
        residue_sequence_number,
        residue_insertion_code,
        x,
        y,
        z,
        occupancy,
        temperature_factor,
        segment_id,
        element_symbol,
        charge;
        pdb_info::PDBInfo{T},
        options) where {T}

    # should we skip this model?
    if ((pdb_info.selected_model != -1) && (pdb_info.selected_model != pdb_info.current_model))
        return
    end

    # hand over parsing to ATOM
    interpret_record(Val(RECORD_TYPE__ATOM), tag, serial_number,
        atom_name, alternate_location_identifier, residue_name,
        chain_id, residue_sequence_number, residue_insertion_code,
        x, y, z, occupancy, temperature_factor, segment_id,
        element_symbol, charge; 
        pdb_info=pdb_info, options=options, is_hetero_atom=true)

    # ensure the residue is marked as a hetero fragement, even if the
    # first atom was a regular one
    set_property!(pdb_info.current_residue, :is_hetero_fragment, true)
    set_property!(pdb_info.current_residue, :is_non_standard, true)

    # recognize water
    if !isnothing(match(r"^OHH|HOH|HHO|H2O|2HO|OH2|SOL|TIP|TIP2|TIP3|TIP4|WAT|D2O$", pdb_info.current_residue.name))
        set_property!(pdb_info.current_residue, :is_water, true)
    end

end

function interpret_record(record_type, tag, record_data...; pdb_info, options)
    push!(pdb_info.records, PDBRecord(tag, record_data))
end

function read_pdb(filename::String, T)
    pdblines = readlines(filename)

    sys = System{T}("")
    mol = Molecule(sys; name="")

    pdb_info = PDBInfo{T}(mol)

    for pl in pdblines
        handle_record(pl, pdb_info)
    end

    sys.name = strip(pdb_info.name)
    mol.name = sys.name

    pdb_info
end

function read_pdb(filename::String)
    read_pdb(filename, Float32)
end