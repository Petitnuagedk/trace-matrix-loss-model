#include "trace-matrix-propagation-model.h"

#include "ns3/log.h"
#include "ns3/node.h"
#include "ns3/mobility-model.h"
#include "ns3/double.h"
#include "ns3/boolean.h"
#include "ns3/enum.h"
#include "ns3/simulator.h"

namespace ns3
{

NS_LOG_COMPONENT_DEFINE("TraceMatrixPropagationLossModel");
NS_OBJECT_ENSURE_REGISTERED(TraceMatrixPropagationLossModel);

TypeId
TraceMatrixPropagationLossModel::GetTypeId()
{
    static TypeId tid =
        TypeId("ns3::TraceMatrixPropagationLossModel")
            .SetParent<PropagationLossModel>()
            .SetGroupName("Propagation")
            .AddConstructor<TraceMatrixPropagationLossModel>()
            .AddAttribute("DefaultLoss",
                          "Loss (dB) used for node pairs with no entry in the trace "
                          "and before InitializeNodes() has been called.",
                          DoubleValue(std::numeric_limits<double>::max()),
                          MakeDoubleAccessor(&TraceMatrixPropagationLossModel::SetDefaultLoss),
                          MakeDoubleChecker<double>())
            .AddAttribute("Symmetric",
                          "If true, assume Loss(a,b) == Loss(b,a) and store only "
                          "the upper triangle. Do not enable this if the underlying "
                          "model captured directional/antenna effects.",
                          BooleanValue(false),
                          MakeBooleanAccessor(&TraceMatrixPropagationLossModel::SetSymmetric),
                          MakeBooleanChecker())
            .AddAttribute("Interpolate",
                          "If true, linearly interpolate loss between the current "
                          "and staged next frame instead of holding the current "
                          "frame constant until the next frame swap.",
                          BooleanValue(false),
                          MakeBooleanAccessor(&TraceMatrixPropagationLossModel::SetInterpolate),
                          MakeBooleanChecker())
            .AddAttribute("SparseThreshold",
                          "In SPARSE storage mode, entries with loss (dB) greater "
                          "than this value are treated as DefaultLoss and not stored. "
                          "Ignored in DENSE mode.",
                          DoubleValue(200.0),
                          MakeDoubleAccessor(&TraceMatrixPropagationLossModel::SetSparseThreshold),
                          MakeDoubleChecker<double>())
            .AddTraceSource("FrameChange",
                            "Fired when a new current frame is committed, with the "
                            "simulation time of the swap.",
                            MakeTraceSourceAccessor(
                                &TraceMatrixPropagationLossModel::m_frameChangeTrace),
                            "ns3::Time::TracedCallback");
    return tid;
}

TraceMatrixPropagationLossModel::TraceMatrixPropagationLossModel()
    : PropagationLossModel(),
      m_dimension(0),
      m_default(std::numeric_limits<double>::max()),
      m_storageMode(StorageMode::SPARSE),
      m_symmetric(false),
      m_interpolate(false),
      m_sparseThreshold(200.0),
      m_currentFrameTime(Seconds(0)),
      m_nextFrameTime(Seconds(0)),
      m_haveNextFrame(false)
{
}

TraceMatrixPropagationLossModel::~TraceMatrixPropagationLossModel()
{
}

void
TraceMatrixPropagationLossModel::InitializeNodes(const std::vector<Ptr<MobilityModel>>& nodes)
{
    m_nodeIdToIndex.clear();
    m_dimension = static_cast<uint32_t>(nodes.size());
    m_lossesDense.clear();
    m_lossesSparse.clear();
    m_nextLossesDense.clear();
    m_nextLossesSparse.clear();
    m_haveNextFrame = false;

    if (m_storageMode == StorageMode::DENSE)
    {
        m_lossesDense.assign(static_cast<size_t>(m_dimension) * m_dimension, m_default);
    }
    // SPARSE mode: leave map empty, absent entries fall back to m_default at lookup time.

    for (uint32_t idx = 0; idx < m_dimension; ++idx)
    {
        Ptr<MobilityModel> mob = nodes[idx];
        NS_ASSERT_MSG(mob, "Null MobilityModel passed to InitializeNodes at index " << idx);
        Ptr<Node> node = mob->GetObject<Node>();
        NS_ASSERT_MSG(node,
                      "MobilityModel at index " << idx
                                                 << " is not aggregated to a Node; "
                                                    "cannot build node-id -> matrix-index map.");
        m_nodeIdToIndex[node->GetId()] = idx;
    }
}

bool
TraceMatrixPropagationLossModel::ValidateFrame(const std::vector<std::vector<double>>& frame) const
{
    if (m_dimension == 0)
    {
        NS_LOG_ERROR("SetCurrentFrame/SetNextFrame called before InitializeNodes(); rejecting.");
        return false;
    }
    if (frame.size() != m_dimension)
    {
        NS_LOG_ERROR("Frame row count (" << frame.size() << ") does not match dimension ("
                                          << m_dimension << "); rejecting frame.");
        return false;
    }
    for (uint32_t i = 0; i < m_dimension; ++i)
    {
        if (frame[i].size() != m_dimension)
        {
            NS_LOG_ERROR("Frame row " << i << " has " << frame[i].size() << " columns, expected "
                                       << m_dimension << "; rejecting frame.");
            return false;
        }
    }
    return true;
}

uint64_t
TraceMatrixPropagationLossModel::PackKey(uint32_t i, uint32_t j)
{
    return (static_cast<uint64_t>(i) << 32) | static_cast<uint64_t>(j);
}

void
TraceMatrixPropagationLossModel::CommitFrame(const std::vector<std::vector<double>>& frame,
                                              std::vector<double>& denseTarget,
                                              std::unordered_map<uint64_t, double>& sparseTarget)
{
    if (m_storageMode == StorageMode::DENSE)
    {
        denseTarget.assign(static_cast<size_t>(m_dimension) * m_dimension, m_default);
        for (uint32_t i = 0; i < m_dimension; ++i)
        {
            uint32_t jStart = m_symmetric ? i : 0;
            for (uint32_t j = jStart; j < m_dimension; ++j)
            {
                double v = frame[i][j];
                denseTarget[i * m_dimension + j] = v;
                if (m_symmetric)
                {
                    denseTarget[j * m_dimension + i] = v;
                }
            }
        }
    }
    else // SPARSE
    {
        sparseTarget.clear();
        for (uint32_t i = 0; i < m_dimension; ++i)
        {
            uint32_t jStart = m_symmetric ? i : 0;
            for (uint32_t j = jStart; j < m_dimension; ++j)
            {
                double v = frame[i][j];
                if (v >= m_sparseThreshold)
                {
                    continue; // treat as default loss, don't store
                }
                sparseTarget[PackKey(i, j)] = v;
                if (m_symmetric && j != i)
                {
                    sparseTarget[PackKey(j, i)] = v;
                }
            }
        }
    }
}

bool
TraceMatrixPropagationLossModel::ValidateSparseEntries(
    const std::vector<SparseLossEntry>& entries) const
{
    if (m_dimension == 0)
    {
        NS_LOG_ERROR("SetCurrentFrameSparse/SetNextFrameSparse called before "
                      "InitializeNodes(); rejecting.");
        return false;
    }
    for (const auto& e : entries)
    {
        if (e.i >= m_dimension || e.j >= m_dimension)
        {
            NS_LOG_ERROR("Sparse frame entry (" << e.i << "," << e.j
                                                 << ") out of range for dimension "
                                                 << m_dimension << "; rejecting frame.");
            return false;
        }
    }
    return true;
}

void
TraceMatrixPropagationLossModel::CommitSparseFrame(
    const std::vector<SparseLossEntry>& entries,
    std::vector<double>& denseTarget,
    std::unordered_map<uint64_t, double>& sparseTarget)
{
    if (m_storageMode == StorageMode::DENSE)
    {
        denseTarget.assign(static_cast<size_t>(m_dimension) * m_dimension, m_default);
        for (const auto& e : entries)
        {
            double v = e.lossDb;
            denseTarget[e.i * m_dimension + e.j] = v;
            if (m_symmetric)
            {
                denseTarget[e.j * m_dimension + e.i] = v;
            }
        }
    }
    else // SPARSE
    {
        sparseTarget.clear();
        for (const auto& e : entries)
        {
            if (e.lossDb >= m_sparseThreshold)
            {
                continue; // treat as default loss, don't store
            }
            double v = e.lossDb;
            sparseTarget[PackKey(e.i, e.j)] = v;
            if (m_symmetric && e.i != e.j)
            {
                sparseTarget[PackKey(e.j, e.i)] = v;
            }
        }
    }
}

bool
TraceMatrixPropagationLossModel::SetCurrentFrameSparse(
    const std::vector<SparseLossEntry>& entries)
{
    if (!ValidateSparseEntries(entries))
    {
        return false;
    }
    CommitSparseFrame(entries, m_lossesDense, m_lossesSparse);
    m_currentFrameTime = Simulator::Now();
    m_frameChangeTrace(m_currentFrameTime);
    return true;
}

bool
TraceMatrixPropagationLossModel::SetNextFrameSparse(const std::vector<SparseLossEntry>& entries,
                                                     Time atTime)
{
    if (!ValidateSparseEntries(entries))
    {
        return false;
    }
    if (atTime <= Simulator::Now())
    {
        NS_LOG_ERROR("SetNextFrameSparse() target time " << atTime
                                                           << " is not in the future; rejecting.");
        return false;
    }
    CommitSparseFrame(entries, m_nextLossesDense, m_nextLossesSparse);
    m_nextFrameTime = atTime;
    m_haveNextFrame = true;
    return true;
}

bool
TraceMatrixPropagationLossModel::SetCurrentFrame(const std::vector<std::vector<double>>& frame)
{
    if (!ValidateFrame(frame))
    {
        return false;
    }
    CommitFrame(frame, m_lossesDense, m_lossesSparse);
    m_currentFrameTime = Simulator::Now();
    m_frameChangeTrace(m_currentFrameTime);
    return true;
}

bool
TraceMatrixPropagationLossModel::SetNextFrame(const std::vector<std::vector<double>>& frame,
                                               Time atTime)
{
    if (!ValidateFrame(frame))
    {
        return false;
    }
    if (atTime <= Simulator::Now())
    {
        NS_LOG_ERROR("SetNextFrame() target time " << atTime
                                                     << " is not in the future; rejecting.");
        return false;
    }
    CommitFrame(frame, m_nextLossesDense, m_nextLossesSparse);
    m_nextFrameTime = atTime;
    m_haveNextFrame = true;
    return true;
}

void
TraceMatrixPropagationLossModel::SetAllLoss(double lossDb)
{
    if (m_dimension == 0)
    {
        NS_LOG_WARN("SetAllLoss() called before InitializeNodes(); no-op.");
        return;
    }
    if (m_storageMode == StorageMode::DENSE)
    {
        std::fill(m_lossesDense.begin(), m_lossesDense.end(), lossDb);
    }
    else
    {
        m_lossesSparse.clear();
        if (lossDb < m_sparseThreshold)
        {
            for (uint32_t i = 0; i < m_dimension; ++i)
            {
                for (uint32_t j = 0; j < m_dimension; ++j)
                {
                    m_lossesSparse[PackKey(i, j)] = lossDb;
                }
            }
        }
    }
}

void
TraceMatrixPropagationLossModel::SetDefaultLoss(double defaultLossDb)
{
    m_default = defaultLossDb;
}

void
TraceMatrixPropagationLossModel::SetStorageMode(StorageMode mode)
{
    m_storageMode = mode;
}

void
TraceMatrixPropagationLossModel::SetSymmetric(bool symmetric)
{
    m_symmetric = symmetric;
}

void
TraceMatrixPropagationLossModel::SetInterpolate(bool interpolate)
{
    m_interpolate = interpolate;
}

void
TraceMatrixPropagationLossModel::SetSparseThreshold(double lossDb)
{
    m_sparseThreshold = lossDb;
}

uint64_t
TraceMatrixPropagationLossModel::GetStoredEntryCount() const
{
    return m_storageMode == StorageMode::DENSE ? m_lossesDense.size() : m_lossesSparse.size();
}

double
TraceMatrixPropagationLossModel::LookupLoss(const std::vector<double>& dense,
                                             const std::unordered_map<uint64_t, double>& sparse,
                                             uint32_t i,
                                             uint32_t j) const
{
    if (m_storageMode == StorageMode::DENSE)
    {
        return dense[i * m_dimension + j];
    }
    auto it = sparse.find(PackKey(i, j));
    return it != sparse.end() ? it->second : m_default;
}

bool
TraceMatrixPropagationLossModel::GetNodeIndex(Ptr<MobilityModel> mob, uint32_t& indexOut) const
{
    if (!mob)
    {
        return false;
    }
    Ptr<Node> node = mob->GetObject<Node>();
    if (!node)
    {
        return false;
    }
    auto it = m_nodeIdToIndex.find(node->GetId());
    if (it == m_nodeIdToIndex.end())
    {
        return false;
    }
    indexOut = it->second;
    return true;
}

double
TraceMatrixPropagationLossModel::DoCalcRxPower(double txPowerDbm,
                                                Ptr<MobilityModel> a,
                                                Ptr<MobilityModel> b) const
{
    if (m_dimension == 0)
    {
        return txPowerDbm - m_default;
    }

    uint32_t ai;
    uint32_t bi;
    if (!GetNodeIndex(a, ai) || !GetNodeIndex(b, bi))
    {
        NS_LOG_WARN("DoCalcRxPower: node pair not found in trace index map; "
                     "falling back to DefaultLoss.");
        return txPowerDbm - m_default;
    }

    double loss = LookupLoss(m_lossesDense, m_lossesSparse, ai, bi);

    if (m_interpolate && m_haveNextFrame)
    {
        Time now = Simulator::Now();
        if (now >= m_currentFrameTime && now <= m_nextFrameTime &&
            m_nextFrameTime > m_currentFrameTime)
        {
            double nextLoss = LookupLoss(m_nextLossesDense, m_nextLossesSparse, ai, bi);
            double span = (m_nextFrameTime - m_currentFrameTime).GetSeconds();
            double elapsed = (now - m_currentFrameTime).GetSeconds();
            double alpha = elapsed / span;
            loss = loss + alpha * (nextLoss - loss);
        }
    }

    return txPowerDbm - loss;
}

int64_t
TraceMatrixPropagationLossModel::DoAssignStreams(int64_t stream)
{
    // This model has no internal RNG use (it replays a deterministic trace),
    // so it consumes no streams. Returning 0 is correct here, not a stub.
    return 0;
}

} // namespace ns3
