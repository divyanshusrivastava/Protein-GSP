function G = get_rig_from_pdb(pdbfilepath, threshold, weighted)
%% Get RIG Model for protein
% For given pdb file path, this function reads the file and return a 
% RIG model for the entered threshold
    pdbstruct = pdbread(pdbfilepath);
    % Get the first model in pdbstruct
    if isfield(pdbstruct, 'Model')
        modelstruct = pdbstruct.Model(1);
    else
        error(message('bioinfo:pdbdistplot:NotAtomicCoordinates'));
    end
    
    coords = zeros(1, 3);
    signal = zeros(1,1);
    seq = '';
    i = 0;
    for atom = modelstruct.Atom
        if strcmp(atom.AtomName,'CA')
            i = i+1;
            seq = strcat(seq, aminolookup(atom.resName));
            coords(i, :) = [atom.X atom.Y atom.Z];
            signal(i) = aahphob_black(aminolookup(atom.resName));
        end
    end

    distance_matrix = squareform(pdist(coords));
    if weighted == true
        weighted_adjacency_matrix = logical(distance_matrix < threshold);
        weighted_adjacency_matrix = weighted_adjacency_matrix - ...
            diag(diag(weighted_adjacency_matrix));
        for i=1:size(weighted_adjacency_matrix,1)
            for j=1:size(weighted_adjacency_matrix,2)
                if weighted_adjacency_matrix(i,j) > 0
                    weighted_adjacency_matrix(i,j) = abs(j-i);
                end
            end
        end
        G = gsp_graph(weighted_adjacency_matrix);
    elseif weighted == false
        unweighted_adjacency_matrix = logical(distance_matrix < threshold);
        unweighted_adjacency_matrix = unweighted_adjacency_matrix - ...
            diag(diag(unweighted_adjacency_matrix));
    
        G = gsp_graph(unweighted_adjacency_matrix);
    end
    
    G = gsp_update_coordinates(G, coords);
    % G = gsp_update_coordinates(G, gsp_compute_coordinates(G, 2));
    
    G.sequence = seq;
    G.signal = signal;
    G = gsp_compute_fourier_basis(G);
    G.pdb_id = pdbstruct.Header.idCode;
end