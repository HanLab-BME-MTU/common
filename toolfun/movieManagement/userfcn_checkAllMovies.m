function userfcn_checkAllMovies(index, value, handles)

if get(handles.checkbox_all, 'Value')
    
    userData = get(handles.figure1, 'UserData');
    for x = setdiff(1:length(userData.MD), userData.id)

        userData.statusM(x).Checked(index) = value;
        set(handles.figure1, 'UserData', userData)

        dfs_checkAllMovies(index, value, handles, x)
    end
end


function dfs_checkAllMovies(index, value, handles, x)

    userData = get(handles.figure1, 'UserData');
    M = userData.dependM;
    
    if value  % If check

            parentI = find(M(index, :));
            parentI = parentI(:)';
            
            if isempty(parentI)

                return
            else
                for i = parentI

                    if userData.statusM(x).Checked(i) || ...
                        (~isempty(userData.package(x).processes_{i}) && ...
                                userData.package(x).processes_{i}.success_ )
                        continue 
                    else
                        userData.statusM(x).Checked(i) = 1;
                        set(handles.figure1, 'UserData', userData)
                        dfs_checkAllMovies(i, value, handles, x)
                    end
                end
            end

    else % If uncheck
            
            subindex = find(M(:,index));
            subindex = subindex(:)';
            
            if isempty(subindex) || ...
                (~isempty(userData.package(x).processes_{index}) ...
                   && userData.package(x).processes_{index}.success_)
                return;
            else
                for i = subindex
                    
                    if userData.statusM(x).Checked(i)
                        
                        userData.statusM(x).Checked(i) = 0;
                        set(handles.figure1, 'UserData', userData)
                        dfs_checkAllMovies(i, value, handles, x)                        
                    end
                end
            end        


    end

