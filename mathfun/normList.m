function [listOfNorms,normedVectors]=normList(vectors)
%calculates the norm of a list of vectors
%
%SYNOPSIS [listOfNorms,normedVectors]=normList(vectors)
%
%INPUT list of vectors (nVectorsXdimension)
%
%OUTPUT listOfNorms: list (nX1) containing the norms of the vectors
%       normedVectors: list (nXdim) containing the normed vectors
%
%c: 1/03 Jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nVectors=size(vectors,1);

listOfNorms=zeros(nVectors,1);
normedVectors=zeros(size(vectors));

listOfNorms=sqrt(vectors(:,1).^2+vectors(:,2).^2+vectors(:,3).^2);
goodVectors=find(listOfNorms);

normedVectors(goodVectors,:)=vectors(goodVectors,:)./(listOfNorms(goodVectors)*[1 1 1]);
