#include "NoduleInfoPerSlice.h"


//////////////////////////////////////////////////////////////////////////
CNoduleInfoPerSlice::CNoduleInfoPerSlice()
{
	m_vNoduleInfo.resize(N_EXPERT_COUNT);
	nSlicePostion = -1;
}


//////////////////////////////////////////////////////////////////////////
CNoduleInfoPerSlice::~CNoduleInfoPerSlice()
{
	m_vNoduleInfo.clear();
}


//////////////////////////////////////////////////////////////////////////
void CNoduleInfoPerSlice::addNoduleInfo(int idxExpert, SNoduleInfo& sNoduleInfo)
{
	if (nSlicePostion >= 0)
	{
		if (sNoduleInfo.nSlicePos != nSlicePostion)
			return;
		else
			m_vNoduleInfo[idxExpert].push_back(sNoduleInfo);
	}
	else
	{
		nSlicePostion = sNoduleInfo.nSlicePos;
		m_vNoduleInfo[idxExpert].push_back(sNoduleInfo);
	}
}

////////////////////////////////////////////////////////////////////////////
//int CNoduleInfoPerSlice::getExpertCount()
//{
//	return (int)m_vNoduleInfo.size();
//}

//////////////////////////////////////////////////////////////////////////
SNoduleInfo* CNoduleInfoPerSlice::getNoduleInfoThisExpert(int idxExpert, int& szNoduleThisSlice)
{
	if (idxExpert >= N_EXPERT_COUNT)
	{
		szNoduleThisSlice = 0;
		return NULL;
	}
	
	szNoduleThisSlice = (int)m_vNoduleInfo[idxExpert].size();
	if (szNoduleThisSlice > 0)
		return &m_vNoduleInfo[idxExpert][0];
	else
		return NULL;
}

//////////////////////////////////////////////////////////////////////////
int CNoduleInfoPerSlice::getSlicePosition()
{
	return nSlicePostion;
}




