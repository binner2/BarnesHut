#pragma once

#include "stdinc.h"
#include <vector>
#include <memory>
#include <cassert>

namespace barnes_hut {

/**
 * @brief High-performance object pool for efficient memory management
 *
 * This pool pre-allocates objects and reuses them to avoid repeated
 * allocation/deallocation overhead. Particularly beneficial for tree nodes
 * in Barnes-Hut algorithm where nodes are constantly created and destroyed.
 *
 * @tparam T Type of objects to pool (must have default constructor and reset() method)
 */
template<typename T>
    requires requires(T obj) {
        { obj.reset() } noexcept -> std::same_as<void>;
    }
class ObjectPool {
public:
    /**
     * @brief Construct object pool with initial capacity
     * @param initial_capacity Number of objects to pre-allocate
     */
    explicit ObjectPool(std::size_t initial_capacity = 1024)
        : next_available_index_(0) {
        reserve(initial_capacity);
    }

    // Delete copy operations (pool manages unique resources)
    ObjectPool(const ObjectPool&) = delete;
    ObjectPool& operator=(const ObjectPool&) = delete;

    // Default move operations
    ObjectPool(ObjectPool&&) noexcept = default;
    ObjectPool& operator=(ObjectPool&&) noexcept = default;

    ~ObjectPool() = default;

    /**
     * @brief Reserve space for additional objects
     * @param capacity Total capacity to reserve
     */
    void reserve(std::size_t capacity) {
        if (capacity > pool_.size()) {
            pool_.reserve(capacity);
        }
    }

    /**
     * @brief Allocate an object from the pool
     * @return Pointer to allocated object (never null)
     *
     * If pool has available objects, returns a reset object.
     * Otherwise, creates a new object and adds it to the pool.
     */
    [[nodiscard]] T* allocate() {
        if (next_available_index_ < pool_.size()) {
            T* obj = pool_[next_available_index_].get();
            obj->reset();
            ++next_available_index_;
            return obj;
        }

        // Need to grow pool
        auto new_object = std::make_unique<T>();
        T* obj_ptr = new_object.get();
        pool_.push_back(std::move(new_object));
        ++next_available_index_;
        return obj_ptr;
    }

    /**
     * @brief Reset pool for reuse (doesn't deallocate memory)
     *
     * This is much faster than clearing and reallocating.
     * Objects remain in memory, ready for next iteration.
     */
    void reset() noexcept {
        next_available_index_ = 0;
    }

    /**
     * @brief Clear pool and free all memory
     */
    void clear() noexcept {
        pool_.clear();
        next_available_index_ = 0;
    }

    /**
     * @brief Shrink pool to fit current usage
     *
     * Reduces memory footprint by removing unused objects.
     * Use this after simulation stabilizes to optimal size.
     */
    void shrink_to_fit() {
        if (next_available_index_ < pool_.size()) {
            pool_.resize(next_available_index_);
            pool_.shrink_to_fit();
        }
    }

    // ========================================================================
    // Statistics and information
    // ========================================================================

    /// Get number of currently allocated objects
    [[nodiscard]] std::size_t allocated_count() const noexcept {
        return next_available_index_;
    }

    /// Get total pool capacity
    [[nodiscard]] std::size_t capacity() const noexcept {
        return pool_.size();
    }

    /// Get number of available objects (without allocation)
    [[nodiscard]] std::size_t available_count() const noexcept {
        return pool_.size() - next_available_index_;
    }

    /// Check if pool is empty
    [[nodiscard]] bool empty() const noexcept {
        return next_available_index_ == 0;
    }

    /// Get memory usage in bytes (approximate)
    [[nodiscard]] std::size_t memory_usage() const noexcept {
        return pool_.size() * sizeof(T) +
               pool_.capacity() * sizeof(std::unique_ptr<T>);
    }

    /**
     * @brief Get pool efficiency (0.0 to 1.0)
     *
     * Higher values mean better memory utilization.
     * Low values suggest the pool is over-allocated.
     */
    [[nodiscard]] double efficiency() const noexcept {
        if (pool_.empty()) return 1.0;
        return static_cast<double>(next_available_index_) /
               static_cast<double>(pool_.size());
    }

private:
    std::vector<std::unique_ptr<T>> pool_;
    std::size_t next_available_index_;
};

// ============================================================================
// Specialized Pool for Tree Nodes with Arena Allocation
// ============================================================================

/**
 * @brief Arena-style allocator for tree nodes
 *
 * Provides even better cache locality by allocating nodes
 * in contiguous memory blocks.
 */
template<typename T>
class ArenaPool {
public:
    static constexpr std::size_t BLOCK_SIZE = 4096; // Cache-friendly block size

    explicit ArenaPool(std::size_t reserve_blocks = 8) {
        blocks_.reserve(reserve_blocks);
    }

    // Delete copy, default move
    ArenaPool(const ArenaPool&) = delete;
    ArenaPool& operator=(const ArenaPool&) = delete;
    ArenaPool(ArenaPool&&) noexcept = default;
    ArenaPool& operator=(ArenaPool&&) noexcept = default;

    ~ArenaPool() = default;

    /**
     * @brief Allocate object from arena
     */
    [[nodiscard]] T* allocate() {
        // Check if current block has space
        if (current_block_index_ < blocks_.size()) {
            auto& block = blocks_[current_block_index_];
            if (next_in_block_ < BLOCK_SIZE) {
                T* obj = &block[next_in_block_];
                obj->reset();
                ++next_in_block_;
                ++total_allocated_;
                return obj;
            }
        }

        // Need new block
        blocks_.emplace_back(BLOCK_SIZE);
        current_block_index_ = blocks_.size() - 1;
        next_in_block_ = 0;
        return allocate();
    }

    /**
     * @brief Reset arena for reuse
     */
    void reset() noexcept {
        current_block_index_ = 0;
        next_in_block_ = 0;
        total_allocated_ = 0;
    }

    /**
     * @brief Clear all memory
     */
    void clear() noexcept {
        blocks_.clear();
        current_block_index_ = 0;
        next_in_block_ = 0;
        total_allocated_ = 0;
    }

    /// Get total allocated objects
    [[nodiscard]] std::size_t allocated_count() const noexcept {
        return total_allocated_;
    }

    /// Get memory usage
    [[nodiscard]] std::size_t memory_usage() const noexcept {
        return blocks_.size() * BLOCK_SIZE * sizeof(T);
    }

private:
    std::vector<std::vector<T>> blocks_;
    std::size_t current_block_index_ = 0;
    std::size_t next_in_block_ = 0;
    std::size_t total_allocated_ = 0;
};

} // namespace barnes_hut
