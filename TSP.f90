
!===============================================================================
!                   SETUP
!===============================================================================

subroutine create_city_array(ncities, lattice_side, seed, &
                            cities_index, cities_coords, distance_matrix)
    !===============================================================================
    ! PURPOSE
    !   Generate a square lattice of side 'lattice_side'. Randomly occupy 
    !   'ncities' sites with a "city", filling 'cities_coords' with their 
    !   coordinates, in order with respect to 'cities_index'.
    !   Fill 'distance_matrix' with the distances between cities; the matrix 
    !   will be symmetric, with null diagonal values.
    !
    ! ARGUMENTS
    !   ncities         (in)    - Number of cities to place on the lattice
    !   lattice_side    (in)    - Side length of the square lattice
    !   seed            (in)    - Random seed (8 integers)
    !   cities_index    (out)   - Python-like index array for cities
    !   cities_coords   (out)   - Coordinates of cities (2, ncities)
    !   distance_matrix (out)   - Distance matrix between cities (symmetric)
    !
    ! NOTES 
    !   - Arrays are stored in column-major order (Fortran default).
    !   - In Python, check that 'lattice_side' is AT LEAST sqrt('ncities')
    !===============================================================================

    implicit none

    ! Local variables
    integer, parameter :: dp=selected_real_kind(13)
    integer :: i, j, x, y, x1, x2, y1, y2, dx, dy
    real(kind=dp) :: rand, dist2
    logical :: occupied

    ! Input arguments
    integer, intent(in) :: ncities
    integer, intent(in) :: lattice_side
    integer, dimension(8), intent(in) :: seed

    ! Output arguments
    integer, dimension(ncities), intent(out)  :: cities_index
    integer, dimension(2, ncities), intent(out) :: cities_coords
    real(kind=dp), dimension(ncities, ncities), intent(out) :: distance_matrix

    !===============================================================================
    ! Main code
    !===============================================================================

    ! Create python-compatible indices for cities
    do i=1,ncities
        cities_index(i) = i-1
    enddo

    ! Choose random coordinates for cities
    call random_seed(put=seed)
    i=1
    do while (i<=ncities)
            ! Generate coordinates
        call random_number(rand)
        x = int(rand*lattice_side)
        call random_number(rand)
        y = int(rand*lattice_side)
            ! Check if a city is already there
        occupied = .false.
        if (i > 1) then
            do j = 1,i-1
                if (cities_coords(1, j) == x .and. cities_coords(2, j) == y) then
                    occupied = .true.
                    exit
                endif
            enddo
        endif
            ! Collect coordinates and update counter if the site is free
        if (.not. occupied) then
            cities_coords(1,i) = x
            cities_coords(2,i) = y
            i = i+1
        endif
    enddo

    ! Create distance matrix (symmetric, null diagonal: compute only upper triangular)
    do i = 1, ncities
        distance_matrix(i, i) = 0.0_dp
        do j = i + 1, ncities
                ! Gather distances
            x1 = cities_coords(1, i)
            x2 = cities_coords(1, j)
            y1 = cities_coords(2, i)
            y2 = cities_coords(2, j)
            dx = x2 - x1
            dy = y2 - y1
            dist2 = real(dx*dx + dy*dy, kind=dp)
                ! Collect distances into matrix
            distance_matrix(i, j) = sqrt(dist2)
            distance_matrix(j, i) = distance_matrix(i, j)
        enddo
    enddo

end subroutine create_city_array

!subroutine asymmetry()
!end subroutine asymmetry

subroutine get_travel_distance(ncities, path, distance_matrix, &
                            travel_distance)
    !===============================================================================
    ! PURPOSE
    !   Calculate the total travel distance of the path resulting from visiting 
    !   the cities located in 'cities_coords', in the order given by 'path'.
    !   The loop must be closed: once arriving in the last city, one must
    !   go back to the starting point.
    !
    ! ARGUMENTS
    !   ncities         (in)    - Number of cities (length of 'path')
    !   path            (in)    - Python-like indices stating the order to travel cities in
    !   distance_matrix (in)    - Matrix with inter-city distances
    !   travel_distance (out)   - Total travel distance
    !
    ! NOTES
    !   Remember to close the loop!
    !===============================================================================

    implicit none

    ! Local variables
    integer, parameter :: dp=selected_real_kind(13)
    integer :: i, from_city, to_city

    ! Input
    integer, intent(in) :: ncities
    integer, dimension(ncities), intent(in) :: path
    real(kind=dp), dimension(ncities,ncities), intent(in) :: distance_matrix

    ! Output
    real(kind=dp), intent(out) :: travel_distance

    !===============================================================================
    ! Main code
    !===============================================================================

    travel_distance = 0.0_dp

    ! Compute distance ("+1" is to convert indexing from python to fortran)
    do i=1,ncities
            ! Select segment extrema
        from_city = path(i)+1
        to_city = path(mod(i, ncities)+1)+1  ! back to start if in the last city
            ! Add segment length
        travel_distance = travel_distance + distance_matrix(from_city,to_city)
    enddo

end subroutine get_travel_distance


!===============================================================================
!                   OPERATORS
!===============================================================================

subroutine swap(ncities, path, distance_matrix, &
                newpath, length_variation)
    !===============================================================================
    ! PURPOSE
    !   Apply the "swap" protocol to a path, swapping the position of two cities.
    !   Return the modified path in 'newpath'.
    !   Compute the "length" variation, through the distance matrix.
    !
    ! ARGUMENTS
    !   ncities             (in)    - Number of cities (length of 'path')
    !   path                (in)    - Python-like indices stating the order to travel cities in
    !   distance_matrix     (in)    - Matrix with inter-city distances
    !   newpath             (out)   - Modified path after swap
    !   length_variation    (out)   - Travel distance variation
    !
    ! NOTES
    !   Two cities are swapped at positions pos1 and pos2
    !===============================================================================
    
    implicit none

    ! Local variables
    integer, parameter :: dp = selected_real_kind(13)
    integer :: pos1, pos2, temp_pos
    integer :: city1, city2, city1_prev, city1_succ, city2_prev, city2_succ
    real(kind=dp) :: rand1, rand2, old_length, new_length

    ! Input
    integer, intent(in) :: ncities
    integer, dimension(ncities), intent(in) :: path
    real(kind=dp), dimension(ncities,ncities), intent(in) :: distance_matrix

    ! Output
    integer, dimension(ncities), intent(out) :: newpath
    real(kind=dp), intent(out) :: length_variation

    !===============================================================================
    ! Main code
    !===============================================================================

    ! Choose two POSITIONS randomly
    call random_number(rand1)
    call random_number(rand2)
        pos1 = int(rand1*ncities)+1
        pos2 = int(rand2*ncities)+1

    ! Make sure they are different
    do while (pos1==pos2)
    call random_number(rand2)
            pos2 = int(rand2*ncities)+1
    end do

    ! Make sure pos1 < pos2
    if (pos1>pos2) then
            temp_pos = pos1
            pos1 = pos2
            pos2 = temp_pos
    end if

    ! Get cities at these positions
    city1 = path(pos1)+1
    city2 = path(pos2)+1

    ! Find adjacent cities
        ! (city1 neighbors)
    if (pos1 == 1) then
            city1_prev = path(ncities)+1
    else
            city1_prev = path(pos1-1)+1
    endif
    if (pos1 == ncities) then
            city1_succ = path(1)+1
    else
            city1_succ = path(pos1+1)+1
    endif
        ! (city2 neighbors)
    if (pos2 == 1) then
            city2_prev = path(ncities)+1
    else
            city2_prev = path(pos2-1)+1
    endif
    if (pos2 == ncities) then
            city2_succ = path(1)+1
    else
            city2_succ = path(pos2+1)+1
    endif

    ! Compute old and new distances
    if (pos2 == pos1 + 1) then  ! Adjacent cities
            old_length = distance_matrix(city1_prev, city1) + &
                        distance_matrix(city1, city2) + &
                        distance_matrix(city2, city2_succ)
            new_length = distance_matrix(city1_prev, city2) + &
                        distance_matrix(city2, city1) + &
                        distance_matrix(city1, city2_succ)
    else   ! Non-adjacent cities
            old_length = distance_matrix(city1_prev, city1) + distance_matrix(city1, city1_succ) + &
                        distance_matrix(city2_prev, city2) + distance_matrix(city2, city2_succ)
            new_length = distance_matrix(city1_prev, city2) + distance_matrix(city2, city1_succ) + &
                        distance_matrix(city2_prev, city1) + distance_matrix(city1, city2_succ)
    endif
        length_variation = new_length - old_length

        ! Create newpath with swapped cities
        newpath = path
        newpath(pos1) = path(pos2)
        newpath(pos2) = path(pos1)

end subroutine swap

subroutine transport(ncities, path, distance_matrix, &
                    newpath, length_variation)
    !===============================================================================
    ! PURPOSE
    !   Apply the "transport" protocol to a path, removing a segment and putting 
    !   it back somewhere else. Return the modified path in 'newpath'.
    !   Compute the "length" variation, through the distance matrix.
    !
    ! ARGUMENTS
    !   ncities             (in)    - Number of cities (length of 'path')
    !   path                (in)    - Python-like indices stating the order to travel cities in
    !   distance_matrix     (in)    - Matrix with inter-city distances
    !   newpath             (out)   - Modified path after transport
    !   length_variation    (out)   - Travel distance variation
    !
    ! NOTES
    !   (none)
    !===============================================================================

    implicit none

    ! Local variables
    integer, parameter :: dp = selected_real_kind(13)
    integer :: start_pos, finish_pos, insert_pos, segment_length
    integer :: city_before_seg, city_after_seg, city_at_insert, city_after_insert
    integer :: first_city_seg, last_city_seg
    integer :: i, new_idx
    real(kind=dp) :: rand, old_length, new_length

    ! Input
    integer, intent(in) :: ncities
    integer, dimension(ncities), intent(in) :: path
    real(kind=dp), dimension(ncities, ncities), intent(in) :: distance_matrix

    ! Output
    integer, dimension(ncities), intent(out) :: newpath
    real(kind=dp), intent(out) :: length_variation

    !===============================================================================
    ! Main code
    !===============================================================================

    ! Choose segment to move and insertion point
    do
        ! Choose segment to move: [start_pos : finish_pos]
        call random_number(rand)
        start_pos = int(rand * (ncities - 1)) + 1
        
        call random_number(rand)
        segment_length = int(rand * min(ncities - start_pos, ncities - 2)) + 1
        finish_pos = start_pos + segment_length - 1
        
        ! Choose insertion point (segment will be inserted before it)
        call random_number(rand)
        insert_pos = int(rand * (ncities - segment_length)) + 1
        
        ! If insert_pos would be in/after the segment, shift it by segment_length
        if (insert_pos >= start_pos) then
            insert_pos = insert_pos + segment_length
        end if
        
        ! Check if this is a non-trivial move
        ! If not, regenerate positions
        if (insert_pos /= start_pos .and. insert_pos /= finish_pos + 1) then
            exit  ! Valid move found, exit loop
        end if
    end do
    
    ! Get relevant cities (Python to Fortran indices)
    first_city_seg = path(start_pos) + 1
    last_city_seg = path(finish_pos) + 1
    
    ! Neighbors of segment in path
    if (start_pos == 1) then
        city_before_seg = path(ncities) + 1
    else
        city_before_seg = path(start_pos - 1) + 1
    end if
    
    if (finish_pos == ncities) then
        city_after_seg = path(1) + 1
    else
        city_after_seg = path(finish_pos + 1) + 1
    end if
    
    ! Neighbors at insertion point
    !   Insert BEFORE insert_pos, so:
    !   - city_at_insert is the city BEFORE insert_pos
    !   - city_after_insert is the city AT insert_pos
    if (insert_pos == 1) then
        city_at_insert = path(ncities) + 1
        city_after_insert = path(1) + 1
    else if (insert_pos > ncities) then  ! if insert_pos = finish_pos + 1 and finish_pos = ncities
        city_at_insert = path(ncities) + 1
        city_after_insert = path(1) + 1
    else
        city_at_insert = path(insert_pos - 1) + 1
        city_after_insert = path(insert_pos) + 1
    end if
    
    ! Old path length
    old_length = distance_matrix(city_before_seg, first_city_seg) + &
                distance_matrix(last_city_seg, city_after_seg) + &
                distance_matrix(city_at_insert, city_after_insert)
    
    ! New path length
    new_length = distance_matrix(city_before_seg, city_after_seg) + &
                distance_matrix(city_at_insert, first_city_seg) + &
                distance_matrix(last_city_seg, city_after_insert)
    
    length_variation = new_length - old_length
    
    ! Build new path indices
    new_idx = 0
    if (insert_pos < start_pos) then
        ! Case 1: Insert before the segment
        ! newpath = [1:insert_pos-1] + [segment] + [insert_pos:start_pos-1] + [finish_pos+1:ncities]
        
        ! Copy [1:insert_pos-1]
        do i = 1, insert_pos - 1
            new_idx = new_idx + 1
            newpath(new_idx) = path(i)
        end do
        
        ! Copy segment [start_pos:finish_pos]
        do i = start_pos, finish_pos
            new_idx = new_idx + 1
            newpath(new_idx) = path(i)
        end do
        
        ! Copy [insert_pos:start_pos-1]
        do i = insert_pos, start_pos - 1
            new_idx = new_idx + 1
            newpath(new_idx) = path(i)
        end do
        
        ! Copy [finish_pos+1:ncities]
        do i = finish_pos + 1, ncities
            new_idx = new_idx + 1
            newpath(new_idx) = path(i)
        end do
        
    else
        ! Case 2: Insert after the segment
        ! newpath = [1:start_pos-1] + [finish_pos+1:insert_pos-1] + [segment] + [insert_pos:ncities]
        
        ! Copy [1:start_pos-1]
        do i = 1, start_pos - 1
            new_idx = new_idx + 1
            newpath(new_idx) = path(i)
        end do
        
        ! Copy [finish_pos+1:insert_pos-1]
        do i = finish_pos + 1, insert_pos - 1
            new_idx = new_idx + 1
            newpath(new_idx) = path(i)
        end do
        
        ! Copy segment [start_pos:finish_pos]
        do i = start_pos, finish_pos
            new_idx = new_idx + 1
            newpath(new_idx) = path(i)
        end do
        
        ! Copy [insert_pos:ncities]
        do i = insert_pos, ncities
            new_idx = new_idx + 1
            newpath(new_idx) = path(i)
        end do
        
    end if
    
end subroutine transport

subroutine reverse(ncities, path, distance_matrix, &
                newpath, length_variation)
    !===============================================================================
    ! PURPOSE
    !   Apply the "reverse" protocol to a path, inverting the order of a segment
    !   in the path. Return the modified path in 'newpath'.
    !   Compute the "length" variation, through the distance matrix.
    !
    ! ARGUMENTS
    !   ncities             (in)    - Number of cities (length of 'path')
    !   path                (in)    - Python-like indices stating the order to travel cities in
    !   distance_matrix     (in)    - Matrix with inter-city distances
    !   newpath             (out)   - Modified path after reversing segment
    !   length_variation    (out)   - Travel distance variation
    !
    ! NOTES
    !   (none)
    !===============================================================================
    
    implicit none
    
    ! Local variables
    integer, parameter :: dp = selected_real_kind(13)
    integer :: pos1, pos2, temp_pos
    integer :: city1, city2, city_before, city_after
    real(kind=dp) :: rand1, rand2, old_length, new_length
    
    ! Input
    integer, intent(in) :: ncities
    integer, dimension(ncities), intent(in) :: path
    real(kind=dp), dimension(ncities, ncities), intent(in) :: distance_matrix
    
    ! Output
    integer, dimension(ncities), intent(out) :: newpath
    real(kind=dp), intent(out) :: length_variation
    
    !===============================================================================
    ! MAIN CODE
    !===============================================================================
    
    ! Select two positions in the path as the segment extrema
        ! Choose two positions randomly
    call random_number(rand1)
    call random_number(rand2)
        pos1 = int(rand1*ncities)+1
        pos2 = int(rand2*ncities)+1
        ! Make sure they are different and not adjacent (it would be a swap protocol in this case)
    do while (pos1==pos2 .or. abs(pos1-pos2)<2)
        call random_number(rand2)
        pos2 = int(rand2*ncities)+1
    end do
        ! Make sure pos1 < pos2
    if (pos1>pos2) then
        temp_pos = pos1
        pos1 = pos2
        pos2 = temp_pos
    end if

    ! Get cities at segment boundaries (Python to Fortran indices)
    city1 = path(pos1)+1      ! First city of segment
    city2 = path(pos2)+1      ! Last city of segment

    ! Find cities adjacent to the segment
        ! City before segment
    if (pos1 == 1) then
        city_before = path(ncities)+1
    else
        city_before = path(pos1-1)+1
    end if
        ! City after segment
    if (pos2 == ncities) then
        city_after = path(1)+1
    else
        city_after = path(pos2+1)+1
    end if

    ! Compute old and new lengths
    old_length = distance_matrix(city_before, city1) + &
                distance_matrix(city2, city_after)
    new_length = distance_matrix(city_before, city2) + &
                distance_matrix(city1, city_after)
    length_variation = new_length - old_length

    ! Build the new path
    newpath(1:pos1-1) = path(1:pos1-1)
    newpath(pos1:pos2) = path(pos2:pos1:-1)
    newpath(pos2+1:ncities) = path(pos2+1:ncities) 

end subroutine reverse

subroutine two_opt(ncities, path, distance_matrix, &
                edge1_start, edge1_end, edge2_start, edge2_end, &
                newpath, length_variation)
    !===========================================================================
    ! PURPOSE
    !   Apply 2-opt to the intersecting cities, removing the intersection.
    !
    ! ARGUMENTS
    !   ncities             (in)  - Number of cities (length of 'path')
    !   path                (in)  - Indices to travel the cities with
    !   distance_matrix     (in)  - Matrix of inter-city distances
    !   edge1_start         (in)  - First city
    !   edge1_end           (in)  - Second city (next to city 1)
    !   edge2_start         (in)  - Third city
    !   edge2_end           (in)  - Fourth city (next to city 3)
    !   newpath             (out) - New instance of the modified path
    !   length_variation    (out) - Resulting variation in travel length
    !
    ! NOTES
    !   City 1 is adjacent in the path to city 2, and the same goes for cities 3 and 4.
    !===========================================================================

    implicit none
    
    ! Local variables
    integer, parameter :: dp = selected_real_kind(13)
    integer :: city_a, city_b, city_c, city_d
    real(kind=dp) :: old_length, new_length
    
    ! Input
    integer, intent(in) :: ncities
    integer, dimension(ncities), intent(in) :: path
    real(kind=dp), dimension(ncities, ncities), intent(in) :: distance_matrix
    integer, intent(in) :: edge1_start, edge1_end, edge2_start, edge2_end

    ! Output variables
    integer, dimension(ncities), intent(out) :: newpath
    real(kind=dp), intent(out) :: length_variation

    !===========================================================================
    ! MAIN CODE
    !===========================================================================
    
    ! Get the cities at edge endpoints (convert to Fortran indices)
    city_a = path(edge1_start) + 1
    city_b = path(edge1_end) + 1
    city_c = path(edge2_start) + 1
    city_d = path(edge2_end) + 1
    
    ! Calculate length variation
    old_length = distance_matrix(city_a, city_b) + distance_matrix(city_c, city_d)
    new_length = distance_matrix(city_a, city_c) + distance_matrix(city_b, city_d)
    length_variation = new_length - old_length
    
    ! Build new path by reversing the segment between edge1_end and edge2_start
    newpath(1:edge1_start) = path(1:edge1_start)
    newpath(edge1_start+1:edge2_start) = path(edge2_start:edge1_end:-1)
    newpath(edge2_start+1:ncities) = path(edge2_end:ncities)
    
end subroutine two_opt

subroutine find_intersecting_edges(ncities, path, cities_coords, &
                                edge1_start, edge1_end, edge2_start, edge2_end, &
                                found_intersection)
    !===========================================================================
    ! PURPOSE
    !   Find two edges in the path that intersect
    !
    ! ARGUMENTS
    !   ncities             (in)  - Number of cities
    !   path                (in)  - Path order (Python indices)
    !   cities_coords       (in)  - City coordinates (2, ncities)
    !   edge1_start         (out) - Position in path where edge 1 starts
    !   edge1_end           (out) - Position in path where edge 1 ends
    !   edge2_start         (out) - Position in path where edge 2 starts
    !   edge2_end           (out) - Position in path where edge 2 ends
    !   found_intersection  (out) - .true. if intersection found
    !===========================================================================

    implicit none
    
    ! Local variables
    integer :: i, j, next_i, next_j
    integer :: x1, y1, x2, y2, x3, y3, x4, y4
    integer :: city_i, city_next_i, city_j, city_next_j

    ! Input
    integer, intent(in) :: ncities
    integer, dimension(ncities), intent(in) :: path
    integer, dimension(2, ncities), intent(in) :: cities_coords
    
    ! Output
    integer, intent(out) :: edge1_start, edge1_end, edge2_start, edge2_end
    logical, intent(out) :: found_intersection
    
    !===========================================================================
    ! MAIN CODE
    !===========================================================================

    ! Initialize outputs
    found_intersection = .false.
    edge1_start = -1
    edge1_end = -1
    edge2_start = -1
    edge2_end = -1
    
    ! Check all pairs of non-adjacent edges
    do i = 1, ncities
        ! Edge i: from path(i) to path(i+1)
        next_i = i + 1
        if (next_i > ncities) next_i = 1
        
        city_i = path(i) + 1      ! Convert to Fortran index
        city_next_i = path(next_i) + 1
        
        x1 = cities_coords(1, city_i)
        y1 = cities_coords(2, city_i)
        x2 = cities_coords(1, city_next_i)
        y2 = cities_coords(2, city_next_i)
        
        ! Only check edges that are at least 2 positions apart
        do j = i + 2, ncities
            next_j = j + 1
            if (next_j > ncities) next_j = 1
            
            ! Don't check if edges share a vertex
            if (next_i == j .or. i == next_j) cycle
            
            city_j = path(j) + 1
            city_next_j = path(next_j) + 1
            
            x3 = cities_coords(1, city_j)
            y3 = cities_coords(2, city_j)
            x4 = cities_coords(1, city_next_j)
            y4 = cities_coords(2, city_next_j)
            
            ! Check if these edges intersect
            call segments_intersect(x1, y1, x2, y2, x3, y3, x4, y4, found_intersection)
            if (found_intersection) then
                edge1_start = i
                edge1_end = next_i
                edge2_start = j
                edge2_end = next_j
                return
            end if
        end do
    end do
    
end subroutine find_intersecting_edges

subroutine segments_intersect(x1, y1, x2, y2, x3, y3, x4, y4, intersect)
    !===========================================================================
    ! PURPOSE
    !   Check if two line segments intersect PROPERLY (not at endpoints,
    !   not when collinear)
    !===========================================================================

    implicit none
    
    ! Local variables
    integer, parameter :: dp = selected_real_kind(13)
    real(kind=dp) :: d1, d2, d3, d4
    real(kind=dp) :: denom, t, u

    ! Input
    integer, intent(in) :: x1, y1, x2, y2, x3, y3, x4, y4

    ! Output
    logical, intent(out) :: intersect
    
    !===========================================================================
    ! MAIN CODE
    !===========================================================================
    
    intersect = .false.
    
    ! Use parametric intersection test
    ! Line 1: P = (x1,y1) + t*(x2-x1, y2-y1)
    ! Line 2: Q = (x3,y3) + u*(x4-x3, y4-y3)
    
    denom = real((x1-x2)*(y3-y4) - (y1-y2)*(x3-x4), kind=dp)
    
    ! If denom = 0, lines are parallel or collinear
    if (abs(denom) < 1.0e-10_dp) then
        intersect = .false.
        return
    end if
    
    ! Calculate parameters t and u
    t = real((x1-x3)*(y3-y4) - (y1-y3)*(x3-x4), kind=dp) / denom
    u = real((x1-x3)*(y1-y2) - (y1-y3)*(x1-x2), kind=dp) / denom
    
    ! Segments intersect if 0 < t < 1 and 0 < u < 1
    if (t > 1.0e-10_dp .and. t < 1.0_dp - 1.0e-10_dp .and. &
        u > 1.0e-10_dp .and. u < 1.0_dp - 1.0e-10_dp) then
        intersect = .true.
    end if
    
end subroutine segments_intersect


!===============================================================================
!                   METROPOLIS
!===============================================================================

subroutine choose_protocol(T, intersected, &
                           protocol)
    !===============================================================================
    ! PURPOSE
    !   Choose the protocol to apply in the current Metropolis step based on  
    !   "thermal" probabilities over the four possible protocols, applying a sort of
    !   discretized inverse transormation method. 2-opt is excluded if there isn't 
    !   at least one intersection.
    !   The thermal probabilities are:
    !       - p(swap)       = T     ||  favoured at high temperatures
    !       - p(transport)  = T     ||  (brute-force)
    !       - p(reverse)    = 1     (always good)
    !       - p(2-opt)      = 1/T   (better at low temperatures, for fine work)
    !
    ! ARGUMENTS
    !   T           (in)  - Equivalent of temperature
    !   intersected (in)  - Boolean; 'true' if there is at least one intersection
    !   protocol    (out) - 0>swap; 1>transport; 2>reverse; 3>2-opt
    !
    ! NOTES
    !   'protocol' will need to be translated into a string once in Python.
    !===============================================================================

    implicit none

    ! Local variables
    integer, parameter :: dp=selected_real_kind(13)
    real(kind=dp) :: rand, c0, c1, c2, c3, p_swap, p_transport, p_reverse, p_2opt, Z

    ! Input
    real(kind=dp), intent(in) :: T
    logical, intent(in) :: intersected

    ! Output
    integer, intent(out) :: protocol

    !===============================================================================
    ! MAIN CODE
    !===============================================================================

    ! CASE: 2-opt feasible
    if (intersected) then

        ! Compute thermal probabilities
        p_swap = T
        p_transport = T
        p_reverse = 1
        p_2opt = 1/T

        ! Normalization factor
        Z = p_swap + p_transport + p_reverse + p_2opt

        ! Compute cumulative "function" (discretized)
        c0 = p_swap/Z
        c1 = c0 + p_transport/Z
        c2 = c1 + p_reverse/Z
        c3 = c2 + p_2opt/Z

        ! Choose randomically
        call random_number(rand)
        if (rand<=c0) then
            protocol = 0
        else if (rand<=c1) then
            protocol = 1
        else if (rand<=c2) then
            protocol = 2
        else
            protocol = 3
        endif

    else   ! CASE: 2-opt non feasible

        ! Compute thermal probabilities
        p_swap = T
        p_transport = T
        p_reverse = 1

        ! Normalization factor
        Z = p_swap + p_transport + p_reverse

        ! Compute cumulative "function" (discretized)
        c0 = p_swap/Z
        c1 = c0 + p_transport/Z
        c2 = c1 + p_reverse/Z

        ! Choose randomically
        call random_number(rand)
        if (rand<=c0) then
            protocol = 0
        else if (rand<=c1) then
            protocol = 1
        else
            protocol = 2
        endif
    endif

end subroutine choose_protocol

subroutine metropolis_step(ncities, cities_coords, path, distance_matrix, temperature, current_travel_length, &
                           newpath, newlength, length_variation, used_protocol, acceptance)
    !===============================================================================
    ! PURPOSE
    !   Compute a single metropolis step for the traveling salesman problem.
    !   Choose a protocol based on the temperature, apply it to the path, and
    !   accept/reject the modification.
    !
    ! ARGUMENTS
    !   ncities                 (in)  - Number of cities (length of 'path')
    !   path                    (in)  - Python-like indices stating the order to travel cities in
    !   distance_matrix         (in)  - Matrix with inter-city distances
    !   T                       (in)  - Equivalent of temperature
    !   current_travel_length   (in)  - Current total travel length
    !   newpath                 (out) - Path after (possible) modification
    !   length_variation        (out) - Variation in the total travel length
    !   used_protocol           (out) - Protocol used in the modification of the path
    !   acceptance              (out) - 1=accepted; 0=refused
    !
    ! NOTES
    !   used_protocol =
    !       - 0     swap
    !       - 1     transport
    !       - 2     reverse
    !       - 3     2-opt
    !===============================================================================

    implicit none

    ! Local variables
    integer, parameter :: dp=selected_real_kind(13)
    integer :: edge1_start, edge1_end, edge2_start, edge2_end
    logical :: found_intersection
    real(kind=dp) :: rand, prob_ratio, proposed_length

    ! Input
    integer, intent(in) :: ncities
    real(kind=dp), intent(in) :: temperature, current_travel_length
    integer, dimension(ncities), intent(in) :: path
    integer, dimension(2,ncities), intent(in) :: cities_coords
    real(kind=dp), dimension(ncities, ncities), intent(in) :: distance_matrix

    ! Output
    integer, intent(out) :: used_protocol, acceptance
    real(kind=dp), intent(out) :: newlength, length_variation
    integer, dimension(ncities), intent(out) :: newpath


    !===============================================================================
    ! MAIN CODE
    !===============================================================================

    ! Choose protocol based on temperature
    call find_intersecting_edges(ncities, path, cities_coords, &
                                edge1_start, edge1_end, edge2_start, edge2_end, &
                                found_intersection)
    call choose_protocol(temperature, found_intersection, used_protocol)

    ! Apply protocol
    if (used_protocol==0) then
        call swap(ncities, path, distance_matrix, newpath, length_variation)
    elseif (used_protocol==1) then
        call transport(ncities, path, distance_matrix, newpath, length_variation)
    elseif (used_protocol==2) then
        call reverse(ncities, path, distance_matrix, newpath, length_variation)
    elseif (used_protocol==3) then
        call two_opt(ncities, path, distance_matrix, edge1_start, edge1_end, edge2_start, edge2_end, newpath, length_variation)
    endif

    ! Recalculate the actual length (adding incremental variation raises roundoff errors and often results in negative final lengths)
    call get_travel_distance(ncities, newpath, distance_matrix, proposed_length)
    length_variation = proposed_length - current_travel_length

    ! Metropolis step
    if (length_variation<=0.0_dp) then  ! (Always accept if better)
        acceptance = 1
        newlength = proposed_length
    else   ! (If worse, evaluate with Metropolis)
        call random_number(rand)
        prob_ratio = exp(-length_variation/temperature)
        if (rand<=prob_ratio) then
            acceptance = 1
            newlength = proposed_length
        else
            acceptance = 0
            newlength = current_travel_length
            newpath = path
            length_variation = 0.0_dp
        endif
    endif

end subroutine metropolis_step

subroutine metropolis_run(nsteps, ncities, cities_coords, initial_path, distance_matrix, temperatures, &
                          paths, lengths, length_variations, &
                          acceptances, total_acceptance_rate, stepwise_acceptance_rate, &
                          used_protocols, protocols_rate, accepted_protocols_rate)
    !===============================================================================
    ! PURPOSE
    !   Compute the full Metropolis run with MC steps done by 'metropolis_step'.
    !   The temperature input is a vector, to allow for different cooling
    !   protocols (that would be chosen and computed in Python). Obtain the vector of
    !   path changes, the length vector, the length variation vector, the acceptance vector, 
    !   the used protocol vector.
    !   ADDED: compute total and step-by-step acceptance rate, proposal rate for each protocol
    !   and usage rate for each protocol.
    !
    ! ARGUMENTS
    !   nsteps                  (in)  - Total Metropolis steps to compute
    !   ncities                 (in)  - Number of cities (length of 'path')
    !   cities_coords           (in)  - Coordinates of the single cities.
    !   initial_path            (in)  - Python-like indices for the starting path
    !   distance_matrix         (in)  - Matrix with inter-city distances
    !   temperatures            (in)  - Vector of temperatures
    !   paths                   (out) - Array of the paths the system followed (dim: ncities x nsteps)
    !   lengths                 (out) - Vector of the path lengths (dim: nsteps)
    !   length_variations       (out) - Vector of the length changes (dim: nsteps)
    !   acceptances             (out) - Vector of the acceptance rates (dim: nsteps)
    !   total_acceptance_rate   (out) - Average acceptance rate for the current run
    !   stepwise_acceptance_rate(out) - Average acceptance rate step-by-ste, in a "cumulative" sense
    !   used_protocols          (out) - Vector of the used protocols
    !   protocols_rate          (out) - How many times each protocol has been proposed, relative to nsteps
    !   accepted_protocols_rate (out) - How many times each protocol has been accepted, relative to accepted steps
    !
    ! NOTES
    !   used_protocol =
    !       - 0     swap
    !       - 1     transport
    !       - 2     reverse
    !       - 3     2-opt
    !   acceptances = 
    !       - 0     refused
    !       - 1     accepted
    !   protocols_rate = [swap_rate, transport_rate, reverse_rate, 2opt_rate]
    !   accepted_protocols_rate = [accepted_swap_rate, accepted_transport_rate, accepted_reverse_rate, accepted_2opt_rate]
    !===============================================================================

    implicit none

    ! Local variables
    integer, parameter :: dp=selected_real_kind(13)
    integer :: i, swap_count, transport_count, reverse_count, opt2_count, accepted_swap_count, accepted_transport_count, accepted_reverse_count, accepted_opt2_count

    ! Input
    integer, intent(in) :: nsteps, ncities
    integer, dimension(ncities), intent(in) :: initial_path
    integer, dimension(2,ncities), intent(in) :: cities_coords
    real(kind=dp), dimension(nsteps), intent(in) :: temperatures
    real(kind=dp), dimension(ncities,ncities), intent(in) :: distance_matrix

    ! Output
    integer, dimension(ncities,nsteps), intent(out) :: paths
    integer, dimension(nsteps), intent(out) :: acceptances, used_protocols
    real(kind=dp), intent(out) :: total_acceptance_rate
    real(kind=dp), dimension(4), intent(out) :: protocols_rate, accepted_protocols_rate
    real(kind=dp), dimension(nsteps), intent(out) :: lengths, length_variations, stepwise_acceptance_rate
    

    !===============================================================================
    ! MAIN CODE
    !===============================================================================

    ! Initialization
    paths(:,1) = initial_path
    length_variations(nsteps) = 0.0_dp
    acceptances(nsteps) = 0
    used_protocols(nsteps) = -1

    call get_travel_distance(ncities, initial_path, distance_matrix, lengths(1))

    ! Run
    do i=1,nsteps-1
        call metropolis_step(ncities, cities_coords, paths(:,i), distance_matrix, temperatures(i), lengths(i), &
                             paths(:,i+1), lengths(i+1), length_variations(i), used_protocols(i), acceptances(i))
    enddo

    ! Other variables: acceptance rates
    total_acceptance_rate = real(sum(acceptances), kind=dp)/nsteps
    stepwise_acceptance_rate(1) = acceptances(1)
    do i=2,nsteps
        stepwise_acceptance_rate(i) = sum(acceptances(1:i))/i
    enddo
    
    ! Other variables: protocols proposals
    
    swap_count = count(used_protocols == 0)
    transport_count = count(used_protocols == 1)
    reverse_count = count(used_protocols == 2)
    opt2_count = count(used_protocols == 3)

    protocols_rate = real([swap_count, transport_count, reverse_count, opt2_count], kind=dp) / nsteps

    ! Other variables: protocols usage
    accepted_swap_count = count(used_protocols == 0 .and. acceptances==1)
    accepted_transport_count = count(used_protocols == 1 .and. acceptances==1)
    accepted_reverse_count = count(used_protocols == 2 .and. acceptances==1)
    accepted_opt2_count = count(used_protocols == 3 .and. acceptances==1)

    accepted_protocols_rate = real([accepted_swap_count, accepted_transport_count, accepted_reverse_count, accepted_opt2_count], kind=dp) / sum(acceptances)

end subroutine metropolis_run

!===============================================================================
!                   SETUP - VERSIONE CON COORDINATE REALI
!===============================================================================

subroutine read_tsplib_coords(ncities, coords_x, coords_y, &
                              cities_index, cities_coords, distance_matrix)
    !===============================================================================
    ! PURPOSE
    !   Read TSPLIB-like files and provide alternative setup tool in the 'create_cit_array' format.
    !   Uses real coordinates and calculates distance matrix.
    !
    ! ARGUMENTS
    !   ncities         (in)    - Number of cities
    !   coords_x        (in)    - X coordinates from TSPLIB
    !   coords_y        (in)    - Y coordinates from TSPLIB
    !   cities_index    (out)   - Python-like index array for cities
    !   cities_coords   (out)   - Coordinates of cities (2, ncities)
    !   distance_matrix (out)   - Distance matrix between cities (symmetric)
    !
    ! NOTES
    !   citiees_coords are now real numbers.
    !===============================================================================

    implicit none

    ! Local variables
    integer, parameter :: dp=selected_real_kind(13)
    integer :: i, j
    real(kind=dp) :: dx, dy, dist2

    ! Input arguments
    integer, intent(in) :: ncities
    real(kind=dp), dimension(ncities), intent(in) :: coords_x, coords_y

    ! Output arguments
    integer, dimension(ncities), intent(out) :: cities_index
    real(kind=dp), dimension(2, ncities), intent(out) :: cities_coords  ! CHANGED TO REAL
    real(kind=dp), dimension(ncities, ncities), intent(out) :: distance_matrix

    !===============================================================================
    ! Main code
    !===============================================================================

    ! Create python-compatible indices for cities
    do i=1,ncities
        cities_index(i) = i-1
    enddo

    ! Store coordinates
    do i=1,ncities
        cities_coords(1,i) = coords_x(i)
        cities_coords(2,i) = coords_y(i)
    enddo

    ! Create distance matrix (symmetric, null diagonal)
    do i = 1, ncities
        distance_matrix(i, i) = 0.0_dp
        do j = i + 1, ncities
            dx = cities_coords(1, j) - cities_coords(1, i)
            dy = cities_coords(2, j) - cities_coords(2, i)
            dist2 = dx*dx + dy*dy
            distance_matrix(i, j) = sqrt(dist2)
            distance_matrix(j, i) = distance_matrix(i, j)
        enddo
    enddo

end subroutine read_tsplib_coords

subroutine find_intersecting_edges_real(ncities, path, cities_coords, &
                                edge1_start, edge1_end, edge2_start, edge2_end, &
                                found_intersection)
    !===========================================================================
    ! PURPOSE
    !   Find two edges in the path that intersect (VERSION WITH REAL COORDS)
    !===========================================================================

    implicit none
    
    ! Local variables
    integer, parameter :: dp=selected_real_kind(13)
    integer :: i, j, next_i, next_j
    real(kind=dp) :: x1, y1, x2, y2, x3, y3, x4, y4  ! CHANGED TO REAL
    integer :: city_i, city_next_i, city_j, city_next_j

    ! Input
    integer, intent(in) :: ncities
    integer, dimension(ncities), intent(in) :: path
    real(kind=dp), dimension(2, ncities), intent(in) :: cities_coords  ! CHANGED TO REAL
    
    ! Output
    integer, intent(out) :: edge1_start, edge1_end, edge2_start, edge2_end
    logical, intent(out) :: found_intersection
    
    !===========================================================================
    ! MAIN CODE
    !===========================================================================

    ! Initialize outputs
    found_intersection = .false.
    edge1_start = -1
    edge1_end = -1
    edge2_start = -1
    edge2_end = -1
    
    ! Check all pairs of non-adjacent edges
    do i = 1, ncities
        next_i = i + 1
        if (next_i > ncities) next_i = 1
        
        city_i = path(i) + 1
        city_next_i = path(next_i) + 1
        
        x1 = cities_coords(1, city_i)
        y1 = cities_coords(2, city_i)
        x2 = cities_coords(1, city_next_i)
        y2 = cities_coords(2, city_next_i)
        
        do j = i + 2, ncities
            next_j = j + 1
            if (next_j > ncities) next_j = 1
            
            if (next_i == j .or. i == next_j) cycle
            
            city_j = path(j) + 1
            city_next_j = path(next_j) + 1
            
            x3 = cities_coords(1, city_j)
            y3 = cities_coords(2, city_j)
            x4 = cities_coords(1, city_next_j)
            y4 = cities_coords(2, city_next_j)
            
            call segments_intersect_real(x1, y1, x2, y2, x3, y3, x4, y4, found_intersection)
            if (found_intersection) then
                edge1_start = i
                edge1_end = next_i
                edge2_start = j
                edge2_end = next_j
                return
            end if
        end do
    end do
    
end subroutine find_intersecting_edges_real

subroutine segments_intersect_real(x1, y1, x2, y2, x3, y3, x4, y4, intersect)
    !===========================================================================
    ! PURPOSE
    !   Check if two line segments intersect (VERSION WITH REAL COORDS)
    !===========================================================================

    implicit none
    
    ! Local variables
    integer, parameter :: dp = selected_real_kind(13)
    real(kind=dp) :: denom, t, u

    ! Input - NOW ALL REAL
    real(kind=dp), intent(in) :: x1, y1, x2, y2, x3, y3, x4, y4

    ! Output
    logical, intent(out) :: intersect
    
    !===========================================================================
    ! MAIN CODE
    !===========================================================================
    
    intersect = .false.
    
    denom = (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4)
    
    if (abs(denom) < 1.0e-10_dp) then
        intersect = .false.
        return
    end if
    
    t = ((x1-x3)*(y3-y4) - (y1-y3)*(x3-x4)) / denom
    u = ((x1-x3)*(y1-y2) - (y1-y3)*(x1-x2)) / denom
    
    if (t > 1.0e-10_dp .and. t < 1.0_dp - 1.0e-10_dp .and. &
        u > 1.0e-10_dp .and. u < 1.0_dp - 1.0e-10_dp) then
        intersect = .true.
    end if
    
end subroutine segments_intersect_real

subroutine metropolis_step_real(ncities, cities_coords, path, distance_matrix, temperature, current_travel_length, &
                           newpath, newlength, length_variation, used_protocol, acceptance)
    !===============================================================================
    ! VERSION WITH REAL COORDINATES
    !===============================================================================

    implicit none

    ! Local variables
    integer, parameter :: dp=selected_real_kind(13)
    integer :: edge1_start, edge1_end, edge2_start, edge2_end
    logical :: found_intersection
    real(kind=dp) :: rand, prob_ratio, proposed_length

    ! Input
    integer, intent(in) :: ncities
    real(kind=dp), intent(in) :: temperature, current_travel_length
    integer, dimension(ncities), intent(in) :: path
    real(kind=dp), dimension(2,ncities), intent(in) :: cities_coords  ! CHANGED TO REAL
    real(kind=dp), dimension(ncities, ncities), intent(in) :: distance_matrix

    ! Output
    integer, intent(out) :: used_protocol, acceptance
    real(kind=dp), intent(out) :: newlength, length_variation
    integer, dimension(ncities), intent(out) :: newpath

    !===============================================================================
    ! MAIN CODE
    !===============================================================================

    ! Choose protocol based on temperature
    call find_intersecting_edges_real(ncities, path, cities_coords, &
                                edge1_start, edge1_end, edge2_start, edge2_end, &
                                found_intersection)
    call choose_protocol(temperature, found_intersection, used_protocol)

    ! Apply protocol (these don't need cities_coords, only distance_matrix)
    if (used_protocol==0) then
        call swap(ncities, path, distance_matrix, newpath, length_variation)
    elseif (used_protocol==1) then
        call transport(ncities, path, distance_matrix, newpath, length_variation)
    elseif (used_protocol==2) then
        call reverse(ncities, path, distance_matrix, newpath, length_variation)
    elseif (used_protocol==3) then
        call two_opt(ncities, path, distance_matrix, edge1_start, edge1_end, edge2_start, edge2_end, newpath, length_variation)
    endif

    ! Recalculate the actual length
    call get_travel_distance(ncities, newpath, distance_matrix, proposed_length)
    length_variation = proposed_length - current_travel_length

    ! Metropolis step
    if (length_variation<=0.0_dp) then
        acceptance = 1
        newlength = proposed_length
    else
        call random_number(rand)
        prob_ratio = exp(-length_variation/temperature)
        if (rand<=prob_ratio) then
            acceptance = 1
            newlength = proposed_length
        else
            acceptance = 0
            newlength = current_travel_length
            newpath = path
            length_variation = 0.0_dp
        endif
    endif

end subroutine metropolis_step_real

subroutine metropolis_run_real(nsteps, ncities, cities_coords, initial_path, distance_matrix, temperatures, seed, &
                          paths, lengths, length_variations, &
                          acceptances, total_acceptance_rate, stepwise_acceptance_rate, &
                          used_protocols, protocols_rate, accepted_protocols_rate)

    !===============================================================================
    ! ('metropolis_run' with real coordinates and seed input
    !===============================================================================

    implicit none

    ! Local variables
    integer, parameter :: dp=selected_real_kind(13)
    integer :: i

    ! Input
    integer, intent(in) :: nsteps, ncities
    integer, dimension(ncities), intent(in) :: initial_path
    integer, dimension(8), intent(in) :: seed
    real(kind=dp), dimension(2,ncities), intent(in) :: cities_coords  ! CHANGED TO REAL
    real(kind=dp), dimension(nsteps), intent(in) :: temperatures
    real(kind=dp), dimension(ncities,ncities), intent(in) :: distance_matrix

    ! Output
    integer, dimension(ncities,nsteps), intent(out) :: paths
    integer, dimension(nsteps), intent(out) :: acceptances, used_protocols
    real(kind=dp), intent(out) :: total_acceptance_rate
    real(kind=dp), dimension(4), intent(out) :: protocols_rate, accepted_protocols_rate
    real(kind=dp), dimension(nsteps), intent(out) :: lengths, length_variations, stepwise_acceptance_rate

    !===============================================================================
    ! MAIN CODE
    !===============================================================================

    ! Initialization
    call random_seed(put=seed)
    paths(:,1) = initial_path
    length_variations(nsteps) = 0.0_dp
    acceptances(nsteps) = 0
    used_protocols(nsteps) = -1

    call get_travel_distance(ncities, initial_path, distance_matrix, lengths(1))

    ! Run
    do i=1,nsteps-1
        call metropolis_step_real(ncities, cities_coords, paths(:,i), distance_matrix, temperatures(i), lengths(i), &
                             paths(:,i+1), lengths(i+1), length_variations(i), used_protocols(i), acceptances(i))
    enddo

    ! Acceptance rates and protocol statistics
    total_acceptance_rate = real(sum(acceptances), kind=dp)/nsteps
    stepwise_acceptance_rate(1) = real(acceptances(1), kind=dp)
    do i=2,nsteps
        stepwise_acceptance_rate(i) = real(sum(acceptances(1:i)), kind=dp)/i
    enddo
    
    protocols_rate(1) = real(count(used_protocols == 0), kind=dp) / nsteps
    protocols_rate(2) = real(count(used_protocols == 1), kind=dp) / nsteps
    protocols_rate(3) = real(count(used_protocols == 2), kind=dp) / nsteps
    protocols_rate(4) = real(count(used_protocols == 3), kind=dp) / nsteps

    accepted_protocols_rate(1) = real(count(used_protocols == 0 .and. acceptances==1), kind=dp) / sum(acceptances)
    accepted_protocols_rate(2) = real(count(used_protocols == 1 .and. acceptances==1), kind=dp) / sum(acceptances)
    accepted_protocols_rate(3) = real(count(used_protocols == 2 .and. acceptances==1), kind=dp) / sum(acceptances)
    accepted_protocols_rate(4) = real(count(used_protocols == 3 .and. acceptances==1), kind=dp) / sum(acceptances)

end subroutine metropolis_run_real