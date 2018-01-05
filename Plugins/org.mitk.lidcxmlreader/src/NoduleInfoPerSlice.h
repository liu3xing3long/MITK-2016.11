#ifndef NoduleInfoPerSlice_h__
#define NoduleInfoPerSlice_h__

#include "Typedefs.h"


class CNoduleInfoPerSlice
{
public:
	CNoduleInfoPerSlice();
	~CNoduleInfoPerSlice();

public:
	void addNoduleInfo(int idxExpert, SNoduleInfo& sNoduleInfo);

public: 
	/* this is depreated, every call getNoduleInfoThisExpert should cover four experts and 
	check the if the return value is empty.
	*/
// 	int				getExpertCount();
	SNoduleInfo*	getNoduleInfoThisExpert(int idxExpert, int& szNoduleThisSlice);
	int				getSlicePosition();

private:
	VVSNoduleInfo m_vNoduleInfo;

private:
	int	nSlicePostion;
};

#endif // NoduleInfoPerSlice_h__
