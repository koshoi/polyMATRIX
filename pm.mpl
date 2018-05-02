with(LinearAlgebra);
with(Student[LinearAlgebra]);

mDegree := proc(mtr :: Matrix) :: Matrix;
	map(degree, mtr);
end proc;

lowestDegree := proc(mtr :: Matrix) :: Vector;
	local mdgr := mDegree(mtr);
	local i := 1, j := 1;
	local icore := 1, jcore := 1;
	local ldegree := mdgr[1, 1];
	for i to RowDimension(mtr) do
		for j to ColumnDimension(mtr) do
			if (mdgr[i, j] < ldegree and mdgr[i, j] >= 0) then
				ldegree := mdgr[i, j];
				icore := i; jcore := j;
			end if;
		end do;
	end do;
	(icore, jcore);
end proc;

isDiagonal := proc(mtr :: Matrix) :: Boolean;
	local i := 1, j := 1;
	local ans := true;
	for i to RowDimension(mtr) do
		for j to ColumnDimension(mtr) do
			if (i <> j) then;
				ans := ans and (mtr[i, j] = 0);
			else;
				ans := ans and (mtr[i, j] <> 0);
			end if;
		end do;
	end do;
	ans;
end proc;

# takes right down corner
simpleSubMatrix := proc(mtr :: Matrix) :: Matrix;
	SubMatrix(mtr, [2..RowDimension(mtr)], [2..ColumnDimension(mtr)]);
end proc;

# first column has only one non-zero element and it is the first one
finishedColumn := proc(mtr :: Matrix) :: Boolean;
	local currCol := Column(mtr, 1);
	local i;
	local ans := true;
	ans := (currCol[1] <> 0);
	for i from 2 to Dimension(currCol) do
		ans := ans and (currCol[i] = 0);
	end do;
	ans;
end proc;

# first row has only one non-zero element and it is the first one
finishedRow := proc(mtr :: Matrix) :: Boolean;
	local currRow := Row(mtr, 1);
	local i;
	local ans := true;
	ans := (currRow[1] <> 0);
	for i from 2 to Dimension(currRow) do
		ans := ans and (currRow[i] = 0);
	end do;
	ans;
end proc;

nonZero := proc(vec :: Vector) :: Int;
	local i;
	for i from 2 to Dimension(vec) do
		if (vec[i] <> 0) then
			return i;
		end if;
	end do;
	return 0;
end proc;

iterationMatrix := proc(amtr :: Matrix) :: Vector(Boolean, Matrix), list;
	local i := 1;
	local j := 1;
	local coef := 1;
	local restrictions := [];
	local iterator;
	local mtr := map(normal, amtr);
	if (not finishedColumn(mtr)) then
		#print("UNFINISHED COLUMN", mtr);
		iterator := Matrix(max(Dimension(mtr)), shape = identity);
		i := nonZero(Column(mtr, 1));
		if (mtr[1,1] <> 0) then
			coef := mtr[1, 1];
			iterator := normal(AddRow(iterator, i, 1, -mtr[i, 1] / coef));
			#print("coef->", coef, "iterator->", iterator);
			return (true, iterator), restrictions;
		else
			iterator := normal(SwapRow(iterator, 1, i));
			return (true, iterator), restrictions;
		end if;
	elif (not finishedRow(mtr)) then
		#print("UNFINISHED ROW", mtr);
		iterator := Matrix(min(Dimension(mtr)), shape = identity);
		j := nonZero(Row(mtr, 1));
		if (mtr[1,1] <> 0) then
			coef := mtr[1, 1];
			iterator := normal(Transpose(AddRow(iterator, j, 1, -mtr[1, j] / coef)));
			#print("coef->", coef, "iterator->", iterator);
			return (false, iterator), restrictions;
		else
			iterator := normal(Transpose(SwapRow(iterator, 1, j)));
			return (false, iterator), restrictions;
		end if;
	else
		if (max(Dimension(mtr)) = 1) then;
			return (true, Matrix(1, 1, 1));
		elif (isDiagonal(mtr)) then;
			return (true, Matrix(max(Dimension(mtr)), shape = identity));
		else
			iterator := iterationMatrix(simpleSubMatrix(mtr));
			return (iterator[1], normal(diagonalConcat(iterator[2])));
		end if;
	end if;
	(true, Matrix(max(Dimension(mtr)), shape = identity));
end proc;

iterate := proc(mtr :: Matrix, lr, iterator ) :: Matrix;
	if (lr) then
		return map(normal, iterator.mtr);
	else
		return map(normal, mtr.iterator);
	end if;
end proc;

diagonalConcat := proc(mtr :: Matrix) :: Matrix;
	local i, j;
	i,j  := Dimension(mtr);
	return <<Matrix(1, 1, 1) | Matrix(1, j)>, <Matrix(i, 1) | mtr>>;
end proc;

swapMatrix := proc(mtr :: Matrix) :: Boolean;
	local ammount := 0;
	local i := 1, j := 1;
	for i to RowDimension(mtr) do
		for j to ColumnDimension(mtr) do
			if (mtr[i, j] <> 0) then;
				ammount := ammount + 1;
			end if;
		end do;
	end do;
	if (ammount = max(Dimension(mtr))) then
		return true;
	else
		return false;
	end if;
end proc;

reverseAction := proc(mtr :: Matrix) :: Matrix;
	local i := 1, j := 1;
	#return mtr;
	if (not swapMatrix(mtr)) then
		for i to RowDimension(mtr) do
			for j to ColumnDimension(mtr) do
				if (i <> j) then;
					mtr[i, j] := mtr[i, j] * (-1);
				end if;
			end do;
		end do;
	end if;
	return mtr;
end proc;

fullCycle := proc(mtr::Matrix) :: list;
	local result := mtr;
	local iterator := (true, Matrix([[1, 1], [1, 1]]), []);
	local right := [];
	local left := [];
	while (not isDiagonal(result) and not isDiagonal(iterator[2])) do
		iterator := iterationMatrix(result);
		result := iterate(result, iterator);
		print("result->", result, iterator);
		if (iterator[1]) then
			left := [op(left), reverseAction(iterator[2])];
		else
			right := [reverseAction(iterator[2]), op(right)];
		end if;
	end do;
	return [op(left), map(normal, result), op(right)];
end proc;

mulMatrixVector := proc(mvec :: list) :: Matrix;
	local result := mvec[1];
	local i;
	for i from 2 to nops(mvec) do
		result := result . mvec[i];
		#print(i, "->", result, mvec[i], result . mvec[i] );
	end do;
	return result;
end proc;
