function projection_length = projectOnPCA(coeff,point_1,point_2)
% Project the vector between points 1 and 2 onto the principal axis
vector_diff = point_1 - point_2;
projection_length = dot(vector_diff, coeff(:,1));
end