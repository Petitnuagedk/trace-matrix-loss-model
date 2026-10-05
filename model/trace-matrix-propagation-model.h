#ifndef TRACE_MATRIX_PROPAGATION_LOSS_MODEL_H
#define TRACE_MATRIX_PROPAGATION_LOSS_MODEL_H

#include "ns3/propagation-loss-model.h"
#include "ns3/traced-callback.h"
#include "ns3/nstime.h"

#include <unordered_map>
#include <vector>
#include <cstdint>

namespace ns3
{

class MobilityModel;
class Node;

/**
 * \brief Replays a pre-captured propagation loss trace as successive matrix
 *        frames instead of recomputing loss from a physical model at every
 *        DoCalcRxPower() call.
 *
 * Design notes / deviations from the naive dense-matrix version:
 *
 *  - Storage can be DENSE or SPARSE (see StorageMode). Real mobile/building
 *    scenarios are dominated by out-of-range pairs sitting at the default
 *    loss, so sparse storage (only entries below a configurable loss
 *    threshold are kept) is usually far cheaper in memory for large N.
 *
 *  - Loss values are stored as double. An earlier version stored loss as
 *    float to halve memory, but that rounding was enough to shift the
 *    computed received signal power by ~1e-6 dB relative to a full
 *    physical-model run with bit-identical RNG streams, which is enough
 *    to flip a marginal Wifi reception success/failure decision. Full
 *    double precision removes that as a source of full-vs-trace
 *    divergence (confirmed via fine-grained per-packet PHY debug logging).
 *
 *  - SYMMETRIC mode stores only the upper triangle when the underlying loss
 *    model has no directional/antenna asymmetry, halving memory again. This
 *    is opt-in, not assumed, because directional antennas or beamforming
 *    break the symmetry assumption.
 *
 *  - Frame validation failures are logged (NS_LOG_ERROR) and rejected
 *    instead of silently discarded. A silently-ignored malformed frame is
 *    worse than a crash: it lets a simulation run to completion on stale
 *    data with no evidence anything went wrong.
 *
 *  - A TracedCallback fires on every accepted frame swap, so frame timing
 *    can be cross-checked against mobility course-change events during
 *    validation instead of relying on ad hoc logging.
 *
 *  - Optional linear interpolation between the current and next frame
 *    (by wall-clock time) avoids the artificial step-function link
 *    up/down transitions that pure zero-order hold introduces at frame
 *    boundaries. This is opt-in via SetInterpolate() because it costs
 *    an extra lookup and requires the next frame to be pre-staged.
 *
 *  - This model is only valid for parameter sweeps that do not feed back
 *    into the loss calculation itself (e.g. TX power). Sweeping antenna
 *    gain, beamforming, or anything angle-of-arrival dependent invalidates
 *    a scalar-per-pair trace; this class does not attempt to detect that
 *    misuse and it is the caller's responsibility.
 */
class TraceMatrixPropagationLossModel : public PropagationLossModel
{
  public:
    enum class StorageMode
    {
        DENSE,
        SPARSE
    };

    /**
     * One (row, column, loss) triplet. Used by the sparse-native frame API
     * so callers building a capture pipeline never have to allocate or
     * populate a dense N x N matrix just to hand it to this class -- the
     * dense vector<vector<double>> API further below still exists for
     * convenience/back-compat, but every entry it doesn't mention still
     * costs a default-fill; the triplet API is the one to use when the
     * capture side already knows most pairs are out of range.
     */
    struct SparseLossEntry
    {
        uint32_t i;
        uint32_t j;
        double lossDb;
    };

    static TypeId GetTypeId();

    TraceMatrixPropagationLossModel();
    ~TraceMatrixPropagationLossModel() override;

    /**
     * Map simulation nodes to matrix row/column indices. Must be called
     * before any SetCurrentFrame()/SetNextFrame() call.
     */
    void InitializeNodes(const std::vector<Ptr<MobilityModel>>& nodes);

    /**
     * Install the frame to use immediately. Returns false and logs an
     * error (does not silently no-op) if the frame is malformed.
     */
    bool SetCurrentFrame(const std::vector<std::vector<double>>& frame);

    /**
     * Stage the frame that will become active at time `atTime`, used only
     * when interpolation is enabled. Returns false on validation failure.
     */
    bool SetNextFrame(const std::vector<std::vector<double>>& frame, Time atTime);

    /**
     * Sparse-native equivalent of SetCurrentFrame(): only entries below the
     * default loss need to be listed; every unlisted pair falls back to
     * DefaultLoss. Unlike the dense API, the caller never builds an N x N
     * structure -- this is the entry point a capture pipeline should use.
     * Indices are validated against the dimension set by InitializeNodes();
     * an out-of-range index rejects the whole frame (logged, not silent).
     */
    bool SetCurrentFrameSparse(const std::vector<SparseLossEntry>& entries);

    /// Sparse-native equivalent of SetNextFrame().
    bool SetNextFrameSparse(const std::vector<SparseLossEntry>& entries, Time atTime);

    void SetAllLoss(double lossDb);
    void SetDefaultLoss(double defaultLossDb);
    void SetStorageMode(StorageMode mode);
    void SetSymmetric(bool symmetric);
    void SetInterpolate(bool interpolate);
    void SetSparseThreshold(double lossDb);

    /// Number of node-pair entries currently stored (diagnostic use).
    uint64_t GetStoredEntryCount() const;

  private:
    double DoCalcRxPower(double txPowerDbm,
                          Ptr<MobilityModel> a,
                          Ptr<MobilityModel> b) const override;
    int64_t DoAssignStreams(int64_t stream) override;

    bool ValidateFrame(const std::vector<std::vector<double>>& frame) const;
    void CommitFrame(const std::vector<std::vector<double>>& frame,
                      std::vector<double>& denseTarget,
                      std::unordered_map<uint64_t, double>& sparseTarget);
    bool ValidateSparseEntries(const std::vector<SparseLossEntry>& entries) const;
    void CommitSparseFrame(const std::vector<SparseLossEntry>& entries,
                            std::vector<double>& denseTarget,
                            std::unordered_map<uint64_t, double>& sparseTarget);
    double LookupLoss(const std::vector<double>& dense,
                       const std::unordered_map<uint64_t, double>& sparse,
                       uint32_t i,
                       uint32_t j) const;
    bool GetNodeIndex(Ptr<MobilityModel> mob, uint32_t& indexOut) const;
    static uint64_t PackKey(uint32_t i, uint32_t j);

    uint32_t m_dimension;
    double m_default;
    StorageMode m_storageMode;
    bool m_symmetric;
    bool m_interpolate;
    double m_sparseThreshold;

    std::unordered_map<uint32_t, uint32_t> m_nodeIdToIndex;

    // Current frame storage (only one of the two is populated, per m_storageMode)
    std::vector<double> m_lossesDense;
    std::unordered_map<uint64_t, double> m_lossesSparse;

    // Next frame storage, used only when m_interpolate is true
    std::vector<double> m_nextLossesDense;
    std::unordered_map<uint64_t, double> m_nextLossesSparse;
    Time m_currentFrameTime;
    Time m_nextFrameTime;
    bool m_haveNextFrame;

    TracedCallback<Time> m_frameChangeTrace;
};

} // namespace ns3

#endif // TRACE_MATRIX_PROPAGATION_LOSS_MODEL_H
