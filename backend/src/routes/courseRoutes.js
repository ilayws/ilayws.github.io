const express = require('express');
const router = express.Router();
const CourseGraph = require('../models/CourseGraph');
const courses = require('../data/courses.json');

// Initialize course graph
const courseGraph = new CourseGraph().initialize(courses);

// Get all courses
router.get('/', (req, res) => {
    try {
        const allCourses = courseGraph.getAllCourses();
        res.json(allCourses);
    } catch (error) {
        res.status(500).json({ error: error.message });
    }
});

// Get course details by ID
router.get('/:id', (req, res) => {
    try {
        const courseDetails = courseGraph.getCourseDetails(req.params.id);
        if (!courseDetails) {
            return res.status(404).json({ error: 'Course not found' });
        }
        res.json(courseDetails);
    } catch (error) {
        res.status(500).json({ error: error.message });
    }
});

// Get prerequisites for a course
router.get('/:id/prerequisites', (req, res) => {
    try {
        const prerequisites = courseGraph.getAllPrerequisites(req.params.id);
        res.json(prerequisites);
    } catch (error) {
        res.status(500).json({ error: error.message });
    }
});

// Get dependent courses
router.get('/:id/dependents', (req, res) => {
    try {
        const dependents = courseGraph.getDependentCourses(req.params.id);
        res.json(dependents);
    } catch (error) {
        res.status(500).json({ error: error.message });
    }
});

// Get available courses based on completed courses
router.post('/available', (req, res) => {
    try {
        const { completedCourses } = req.body;
        if (!Array.isArray(completedCourses)) {
            return res.status(400).json({ error: 'completedCourses must be an array' });
        }
        const available = courseGraph.getAvailableCourses(new Set(completedCourses));
        res.json(available);
    } catch (error) {
        res.status(500).json({ error: error.message });
    }
});

// Find path to target course
router.post('/path', (req, res) => {
    try {
        const { targetCourse, completedCourses } = req.body;
        if (!targetCourse || !Array.isArray(completedCourses)) {
            return res.status(400).json({ error: 'Invalid request body' });
        }
        const path = courseGraph.findShortestPath(targetCourse, new Set(completedCourses));
        res.json(path);
    } catch (error) {
        res.status(500).json({ error: error.message });
    }
});

// Filter courses
router.post('/filter', (req, res) => {
    try {
        const criteria = req.body;
        const filteredCourses = courseGraph.filterCourses(criteria);
        res.json(filteredCourses);
    } catch (error) {
        res.status(500).json({ error: error.message });
    }
});

module.exports = router; 