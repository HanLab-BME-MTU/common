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
[nVectors,nDims]=size(vectors);

listOfNorms=zeros(nVectors,1);
normedVectors=zeros(nVectors,nDims);

listOfNorms=sqrt(sum(vectors.^2,2));
goodVectors=find(listOfNorms);

normedVectors(goodVectors,:)=vectors(goodVectors,:)./(repmat(listOfNorms(goodVectors),[1,nDims]));
