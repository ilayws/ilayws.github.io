const express = require('express');
const cors = require('cors');
const helmet = require('helmet');
const morgan = require('morgan');
const coursesRouter = require('./routes/courses');

const app = express();
const port = process.env.PORT || 3000;

// Middleware
app.use(helmet()); // Security middleware
app.use(cors()); // Enable CORS for all routes
app.use(morgan('dev')); // Request logging
app.use(express.json()); // Parse JSON request bodies

// Routes
app.use('/api/courses', coursesRouter);

// Error handling middleware
app.use((err, req, res, next) => {
    console.error(err.stack);
    res.status(500).json({
        error: 'Internal Server Error',
        message: process.env.NODE_ENV === 'development' ? err.message : undefined
    });
});

// 404 handler
app.use((req, res) => {
    res.status(404).json({
        error: 'Not Found',
        message: 'The requested resource was not found'
    });
});

// Start server
app.listen(port, () => {
    console.log(`Server is running on port ${port}`);
}); 